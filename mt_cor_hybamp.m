%MT_COR_HYBAMP - Calculates hybridization and amplification model correction factors
%
%  PROBES = MT_COR_HYBAMP(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure for which to calculate cor factors
%   VARARGIN        'use_amplification' : Enable amplification correction
%                   'm_estimation' : Use m-estimator to estimate model params
%                   'simple'    :    Use simplified sequence model
%                   'reduce_probes' : Fit model with subset of probes. Should be followed by
%                                     probe ratio (between 0 and 1) that should be used.
%
% OUTPUT
%   PROBES		    Probes structure with cor. factors added
%
% DESCRIPTION
% Calculates hybridization and amplification correction factors (after application of ealier 
% corrections) and adds them to the probes structure as fields 'seq_correction' and 'amp_correction'
%
% SEE ALSO
% MT_BG_EST, MT_COR_QQ, MT_COR_IMAGE, MT_REAL_SIGNAL
%
% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function probes = mt_cor_hybamp(probes,varargin)
nprobe = size(probes.pm,2);
narray = size(probes.pm,1);
ngene = length(probes.name);

range = 1:narray;

for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'm_estimation',
            m_estimation = 1;
         case 'simple',
            simple = 1;
         case 'use_amplification',
            amplif = 1;
         case 'ampli_splines',
            ampli_splines = 1;
         case 'range'
            i = i + 1;
            range = varargin{i};
         case 'reduce_probes',
            i = i + 1;
            reduce_probes = varargin{i};
      end;
   end;
end;

signal = 2.^mt_real_signal(probes,'no_qq','no_sequence');
min_signal = min(signal);
bg = min([probes.seqbg(:)'; (min_signal(:)'-1)])';
o = probes.seqbg_factors(1);
bg = bg - o;
med_sig = median(signal)';

x = med_sig -bg - o;

%prevent negative x
r = x < 0;
bg(r) =bg(r) + x(r);
x(r) = 0;

%prevent negative b
r = bg < 0;
bg(r) = 0;  

if(exist('simple'))
   C = mt_create_seq_model(probes,'basis',1);
else
   C = mt_create_seq_model(probes,'basis',2);
end;
nvars = size(C,2)

if(exist('amplif'))
   if(exist('ampli_splines'))
      D = mt_create_amp_model(probes,'splines');
   else
      D = mt_create_amp_model(probes);
   end;      
   nvars_amp = size(D,2)
else
   D = [];
end;

opts = [1e-3 1e-4 1e-6 100 1e-14];

if(exist('reduce_probes'))
    fil = randperm(size(C,1));
    fil = fil(1:ceil(reduce_probes * size(C,1)));
    C = C(fil,:);
    if(~isempty(D))
        D = D(fil,:);
    end;
    signal = signal(:,fil);
    bg = bg(fil);
    x = x(fil);
    nprobe = length(fil);
end;



for ri = 1:length(range)
   fprintf(1,'Array %d: ',ri)
   i = range(ri);
   pms = log2(signal(i,:))';
   weights = ones(nprobe,1, 'single');
   pre_res = zeros(nprobe,1, 'single');
   converged = 0;

   if(exist('amplif'))
      xk = [o; zeros(nvars,1, 'single'); zeros(nvars,1, 'single'); zeros(nvars_amp,1, 'single')];
   else
      xk = [o; zeros(nvars,1, 'single'); zeros(nvars,1, 'single')];
   end;
   while(~converged)

      F = @(vars)calc_score(vars,pms,C,weights,bg,x,i,D);
      [xk,info] = marquardt(F,xk,opts);
      if(info(6) <= 0)
        i
        info(6)
        info
      end;
      fprintf(1,'*');
      residual = F(xk);
      residual = residual ./ weights;

      tol_val = mean(abs(residual - pre_res));
      if(~exist('m_estimation') || (tol_val < 0.01))
         converged = 1;
      else
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
      if(~exist('m_estimation'))
         converged = 1;
      end;
   end;
   probes.seq_factors(i,:) = xk';
   probes.seq_norm(i) = sum(residual.^2);
   fprintf(1,'\n');
end;
%probes.seq_factors(narray+1,:) = 0;

if(exist('reduce_probes'))
    clear signal C D 
    signal = 2.^mt_real_signal(probes,'no_qq','no_sequence');
    min_signal = min(signal);
    bg = min([probes.seqbg(:)'; (min_signal(:)'-1)])';
    o = probes.seqbg_factors(1);
    bg = bg - o;
    med_sig = median(signal)';

    x = med_sig -bg - o;

    %prevent negative x
    r = x < 0;
    bg(r) =bg(r) + x(r);
    x(r) = 0;

    %prevent negative b
    r = bg < 0;
    bg(r) = 0;  

    if(exist('simple'))
        C = mt_create_seq_model(probes,'basis',1);
    else
        C = mt_create_seq_model(probes,'basis',2);
    end;
    nvars = size(C,2)

    if(exist('amplif'))
        if(exist('ampli_splines'))
            D = mt_create_amp_model(probes,'splines');
        else
            D = mt_create_amp_model(probes);
        end;      
        nvars_amp = size(D,2)
    else
        D = [];
    end;
end;

sc = calc_correction(probes.seq_factors,C,o,bg,x,signal,D);
probes.ampseq_correction = log2(signal) - log2(signal - sc);
%if(exist('amplif'))
%   sc2 = calc_correction(probes.seq_factors,C,o,bg,x,signal,D);
%   probes.amp_correction = log2(signal - sc) - log2(signal -sc2);
%end;
return;

function c = calc_correction(seq_factors,C,o,bg,x,signal,D)
[nprobe,nvars] = size(C);
nvars_amp = size(D,2);

narray = size(seq_factors,1);

est_os = seq_factors(:,1);
est_bgseqf = 2.^(C * seq_factors(:,2:(nvars + 1))');
est_fgseqf = 2.^(C * seq_factors(:,(2+nvars):(1 + 2 * nvars))');
if(~isempty(D))
   est_fgseqf = est_fgseqf .* (2.^(D * seq_factors(:,(2 + 2 * nvars):(1 + 2 * nvars + nvars_amp))'));
end;

est_b = repmat(bg,1,narray)' .* est_bgseqf';
est_o = repmat(est_os,1,nprobe);
est_x_w = (signal - est_b - est_o);

%prevent negative x
r = est_x_w < 0;
est_b(r) = est_b(r) + est_x_w(r);
est_x_w(est_x_w < 0) = 0;

%prevent negative b
r = est_b < 0;
est_b(r) = 0;

est_x = est_x_w ./ est_fgseqf';

c = (est_o - o) + (est_b - repmat(bg,1,narray)') + (est_x_w - est_x); 
return;


function [f,J] = calc_score(vars,pms,C,weights,bg,x,arraynr,D)
   [nprobe,nvars] = size(C);
   nvars_amp = size(D,2);
   optic = vars(1);
   bgseqf = 2.^(C * vars(2:(nvars + 1)));
   fgseqf = 2.^(C * vars((2 + nvars):(1 + 2 * nvars)));
   if(~isempty(D))
      fgseqf = fgseqf .* 2.^(D * vars((2 + 2 * nvars):(1 + 2 * nvars + nvars_amp)));
   end;
   
   %prevent infinite
   bgseqf(bgseqf > 2.^500) = 2.^500;
   fgseqf(fgseqf > 2.^500) = 2.^500;
   l = optic + bg.*bgseqf + x .* fgseqf;
   
   %prevent too low values
   l(l <= 0.00001) = 0.00001;

   f = weights .* (pms - log2(l));
   fprintf(1,'.')
   if(nargout > 1)
      k = - weights ./ (l .* log(2));
      Jb = - weights .* bg .* bgseqf ./ l;
      Jf = - weights .* x .* fgseqf ./ l;
      if(~isempty(D))
         J = [k C C D];
         for i = 2:(nvars + 1)
            J(:,i) = J(:,i) .* Jb;
         end;
         for i = (nvars + 2):(nvars * 2 + nvars_amp + 1)
            J(:,i) = J(:,i) .* Jf;
         end;
         %J = [k (C .* (Jb * ones(1,nvars))) (C .* (Jf * ones(1,nvars))) (D .* (Jf * ones(1,nvars_amp)))];
      else

         J = [k C C];
         for i = 2:(nvars + 1)
            J(:,i) = J(:,i) .* Jb;
         end;
         for i = (nvars + 2):(nvars * 2 + 1)
            J(:,i) = J(:,i) .* Jf;
         end;
         %J = [k (C .* (Jb * ones(1,nvars))) (C .* (Jf * ones(1,nvars)))];
      end;
   end;
return;



