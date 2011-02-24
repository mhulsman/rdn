%MT_BG_EST - Estimates the optical and background signal for each probe
%
%  PROBES = MT_BG_EST(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure for which to calculate background estimate.
%   VARARGIN        'reg_max',x   : how many times a pos. residue for (bgest - minsig) 
%                                   should be weighted over a neg. residue (default: 10)
%                   'm_estimation': Use Huber M-estimator
%                   'reduce_probes': Reduce number of probes used in model fitting to increase speed
%                                    and reduce memory usage. Should be followed
%                                    by factor between 0 and 1 to indicate relative number
%                                    of probes to use (1 = all probes, 0 = no probes).
%
% OUTPUT
%   PROBES		    Probes structure with bg estimate fields added
%
% DESCRIPTION
% Estimates the optical and background signal for each probe using the minimum signal over
% all arrays, and fitting a probe sequence model
%
% SEE ALSO
% MT_COR_QQ, MT_COR_IMAGE, MT_COR_HYBAMP, MT_REAL_SIGNAL

% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function probes = mt_bg_est(probes,varargin)

nprobe = size(probes.pm,2);
narray = size(probes.pm,1);
ngene = length(probes.name);

reg_max = 10;

for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'm_estimation',
            m_estimation = 1;
         case 'reg_max',
            i = i + 1;
            reg_max = varargin{i};
         case 'reduce_probes',
            i = i + 1;
            reduce_probes = varargin{i};
      end;
   end;
end;

%sequence model (nprobe * nvar)
C = mt_create_seq_model(probes,'basis',2);
nvars = size(C,2);

%min signal used to estimate bg from sequence
signal = mt_real_signal(probes);
if(exist('reduce_probes'))
    fil = randperm(size(C,1));
    fil = fil(1:ceil(reduce_probes * size(C,1)));
    C = C(fil,:);
    signal = signal(:,fil);
    nprobe = length(fil);
end;


pms = min(signal)';

% m-estimator weights
weights = ones(nprobe,1,'single');

% residual in last iteration
pre_res = zeros(nprobe,1,'single');

% m-estimator reiterative least squares converged indication 
converged = 0;

% options for marquardt least square optimizer
opts = [1e-3 1e-4 1e-8 100 1e-14];
%initial conditions (optical bg, hyb. model params)
xk = [20; zeros(nvars,1,'single')];
fprintf(1,'Estimate background');
while(~converged)
   fprintf(1,'*'); 
   
   %perform optimization and calculate residuals
   F = @(vars)calc_score(vars,pms,C,weights,reg_max);
   [xk] = marquardt(F,xk,opts);
   residual = F(xk);
  
   %get original unweighted residuals back
   residual = residual ./ weights;

   %unsymmetric least square weighting correction
   neg_idx = find(residual < 0);
   residual(neg_idx) = residual(neg_idx) ./ reg_max;

   %difference with residual in last iteration
   tol_val = mean(abs(residual - pre_res));
   
   %check for convergence
   if(tol_val < 0.01 || ~exist('m_estimation'))
      converged = 1;
   else
      %calculate new m-estimator weights
      scale = median(abs(residual))/0.6745;
      res_scale = residual / scale;
      select = (abs(res_scale) <= 1.345);
      if(sum(select) == length(res_scale))
         converged = 1;
      else
         weights(select) = 1;
         weights(~select) = 1.345./abs(res_scale(~select));
      end;
      pre_res = residual;
   end;
end;
fprintf(1,'\n');

%save model parameters
probes.seqbg_factors = xk';

clear signal
%save actual background calculated using model
if(exist('reduce_probes'))
    clear C
    C = mt_create_seq_model(probes,'basis',2);
end;
probes.seqbg = calc_bg(probes.seqbg_factors',C);

return;

function [f,J] = calc_score(vars,pms,C,weights,reg_max)
   %optical background
   optic = vars(1);
   
   %background hybridization model values
   bgseqf = 2.^(C * vars(2:end));
   
   %protect against infinity problems
   bgseqf(bgseqf > 2.^50) = 2.^50;

   %calculate optical + hybridization background
   l = optic + bgseqf;

   %protect against log2 NaN problems
   l(l < 0.00001) = 0.00001;

   % calculate assymetric weighted residue
   f = weights .* (pms - log2(l));
   neg_idx = uint32(find(f < 0));
   f(neg_idx) = reg_max * f(neg_idx);
   fprintf(1,'.');
   %calculate Jacobian if necessary
   if(nargout > 1)
      %combine weights and assymetric residue weighting
      nweights = weights;
      nweights(neg_idx) = reg_max * weights(neg_idx);

      [nprobe,nvars] = size(C);

      %adapt weighted Cw matrix with assymetric residue weighting
      for j = 1:size(C,2)
         C(:,j) = nweights .* C(:,j);
      end;

      %optical parameter Jacobian values
      k = nweights ./ (l .* log(2));
      
      %hybridization parameters Jacobian value part
      k2 = bgseqf ./ l;
      J = [-k -C];
      for r = 1:nvars
        J(:,r+1) = J(:,r+1) .* k2;
      end;
   end;
return;


function bg = calc_bg(vars,C)
   bgseqf = 2.^(C * vars(2:end));
   %sum optical background and background seq. model values
   bg = vars(1) + bgseqf;
return;
