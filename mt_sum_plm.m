%MT_COR_HYBAMP - Summarizes (corrected) probe signal over probesets
%
%  E = MT_SUM_PLM(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure for which to summarize
%   VARARGIN        'm_estimation' : Use m-estimator to estimate summarization parameters
%                   'no_bgadapt'    :   Do not backscale array factors using average background
%                   'no_remove'     :   Do not remove probes which perform bad
%                   'keep_probe_info':  Return an augmented probes structure instead of an expression structure
%
% OUTPUT
%   E		    Expression structure with summarization fields added
%
% DESCRIPTION
% Summarizes the corrected probe signal over probesets.
%
% SEE ALSO
% MT_BG_EST, MT_COR_QQ, MT_COR_IMAGE, MT_COR_HYBAMP
%
% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function r = mt_sum_plm(probes,varargin)

for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'use_batch',
            use_batch = 1;
         case 'm_estimation',
            m_estimation = 1;
         case 'no_bgadapt',
            no_bgadapt = 1;
         case 'no_remove', 
            no_remove = 1;
         case 'keep_probe_info',
            keep_probe_info = 1;
      end;
   end;
end;


%common variables
nprobe = size(probes.pm,2);
narray = size(probes.pm,1);
ngene = length(probes.name);


%handle possible batch labels (add batch factor)
if(isfield(probes,'lab'))
   labels = unique(probes.lab);
   nlab = length(labels);
   lab_idx = {};
   for i = 1:nlab
      lab_idx{i} = find(probes.lab == labels(i));
   end;
end;

%determine number of probes per gene
indexm = sparse(double(probes.ind),1:nprobe,ones(nprobe,1));
nprobe_gene = sum(indexm,2);
nprobe_gene_u = unique(nprobe_gene);

mnprobe_gene = full(max(nprobe_gene))
%precalculate multiplication matrices
%opt_matrices = {};
fprintf(1,' Preparing optimization matrices');
for i = 1:mnprobe_gene
   fprintf(1,'.');
   nx = i;
   
   x = zeros(nx + narray,1);
   y = zeros(nx * narray,1);

   %probe effects
   x(1:(nx * narray)) = reshape(repmat(1:nx,narray,1),[nx * narray,1]);
   y(1:(nx * narray)) = 1:(nx * narray);
   start = nx * narray;
   
   %array effect
   x(start + (1:(nx *narray))) = nx + reshape(repmat(1:narray,1,nx),[nx * narray,1]);
   y(start + (1:(nx *narray))) = 1:(nx * narray);

   s = ones(2 * nx * narray,1);
   opt_matrices{nx} = sparse(y,x,s);
end;
fprintf(1,'\n')

signal = mt_real_signal(probes);


if(exist('keep_probe_info'))
    r = probes;
else
    r = mt_empty_e(probes);
end;

%ASSUMES THAT PROBES ARE ORDERED ON PROBE NR!
if(sum(~(probes.ind == sortrows(probes.ind,1))) > 0)
   fprintf(1,'ERROR: probes hould be ordered per gene and probe-nr');
else
   fprintf(1,' Calculating PLMs');
   r.probe_factors = zeros(1,nprobe, 'single');
   r.removed_factors = zeros(1,nprobe, 'single');
   r.array_factors = zeros(narray,ngene, 'single');
   r.overall_factors = zeros(1,ngene, 'single');
   r.probeset_bg = zeros(1,ngene, 'single');
   if(exist('use_batch') && isfield(probes,'lab'))
      r.batch_factors = zeros(nlab,ngene, 'single');
   end;

   pos = 1;
   
   opts = [1e-3 1e-4 1e-8 100 1e-14];

   min_signal = 2.^min(signal);
   sbg = min([probes.seqbg(:)'; (min_signal(:)'-1)])';

   for i = 1:ngene
      if ((mod(i,1000)==0)), fprintf(1,'.'); end;
      nx = nprobe_gene(i);
      tnx = nx;
      d = signal(:,pos:(pos + nx - 1));
      bg = repmat(sbg(pos:(pos + nx -1))',narray,1);

      converged = 0;
      
      weights = ones(size(bg), 'single');
      pre_res = zeros(size(bg(:)),'single');
      active_probes = (1:nx)';
      removed_probes = [];
      loops = 0;
      
      x = [mean(d)'; zeros(narray,1)];
      while(~converged)
         loops = loops + 1;
         F = @(vars)calc_score(vars,bg,weights,d,single(full(opt_matrices{tnx})));
         [x,info] = marquardt(F,x,opts);
         if(info(6) <= 0)
            i
            info(6)
            info
         end;
         residual = F(x);

         q = find(weights(:) > 0);
         residual(q) = residual(q) ./ weights(q);
         tol_val = mean(abs(residual - pre_res));
         if(~exist('m_estimation') || tol_val < 0.01)
            converged = 1;
            %loops
         else
            %huber weights estimation
            oweights = weights;
            scale = median(abs(residual))/0.6745;
            res_scale = residual / scale;
            select = (abs(res_scale) <= 1.345);
            if(sum(select) == length(res_scale))
               converged = 1;
            else
               weights(select) = 1;
               weights(~select) = 1.345./abs(res_scale(~select));
            end;

            %remove probes with more than 33% weights < 0.9
            tfs = find(sum(weights < 0.9) >= (narray / 3));
            if(~(tnx - length(tfs) < 5) && ~exist('no_remove'))
               removed_probes = [removed_probes; active_probes(tfs)];
               active_probes(tfs) = [];
               weights(:,tfs) = [];
               bg(:,tfs) = [];
               x(tfs) = [];
               d(:,tfs) = [];
               residual = reshape(residual,[narray tnx]);
               residual(:,tfs) = [];
               residual = residual(:);
               tnx = tnx - length(tfs);
            end;
            pre_res = residual;
         end;
      end;

      probef= x(1:tnx)';
     
      arrayf = x(tnx + (1:narray));
      overallf = mean(arrayf) + mean(probef);
      probef = probef - mean(probef);
      arrayf = arrayf - mean(arrayf);
      
      if(~exist('no_bgadapt'))
         ref = log2(2.^overallf + mean(bg(:)));
         arrayf = log2(2.^(overallf + arrayf) + mean(bg(:))) - ref;
      end;
      
      if(isfield(probes,'lab') && exist('use_batch'))
         for j = 1:nlab
            batchf(j) = median(arrayf(lab_idx{j}));
            arrayf(lab_idx{j}) = arrayf(lab_idx{j}) - batchf(j);
         end;
         r.batch_factors(:,i) = batchf;
      end; 
    
      probef2 = [];
      probef2(active_probes) = probef;
      probef2(removed_probes) = -1;
      r.probe_factors(pos:(pos + nx - 1)) = probef2;
      r.overall_factors(i) = overallf;
      r.array_factors(:,i) = arrayf;
      r.removed_factors(pos:(pos + nx -1)) = (probef2 == -1);
      r.probeset_bg(i) = mean(bg(:));
      pos = pos + nx;
   end;
   fprintf(1,'\n');
end;



function [f,J] = calc_score(vars,bg,weights,pms,optmat)
   [narray,nprobe] = size(pms);
   pa = 2.^vars(1:nprobe);
   af = 2.^vars(nprobe + (1:narray));
   x = af * pa';
   l = bg + x;
   f = weights .* (pms - log2(l));
   f = f(:);
   
   k = -x .* weights ./ l;
   if(nargout > 1)
       J = optmat .* (k(:) * ones(1, nprobe + narray));
   end;
return

