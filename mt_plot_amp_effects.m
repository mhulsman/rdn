%MT_PLOT_SEQ_EFFECTS - Plots amplification effects for all arrays
%
%  MT_PLOT_AMP_EFFECTS(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure for which to calculate cor factors
%   VARARGIN        'steps'     : Number of steps along probeset sequence. Should equal
%                                 default number of probes. Default: 11
%                   'difference'  : Show difference values w.r.t mean   
%                                     
%
% DESCRIPTION
%   Show amplification effect differences between arrays, by taking
%   median expression for first probes in probesets, second probes in probests,etc.
%
% SEE ALSO
% MT_COR_HYBAMP, MT_PLOT_ARRAY, MT_PLOT_SEQ_EFFECTS, MT_PLOT_GENE
%
% (c) Marc Hulsman, 2009
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands


function mt_plot_amp_effects(probes,varargin)

   % Default parameters
   steps = 11;
   %should be enough probes in a probeset
   ncount = mt_countprobes(probes);
   probes = mt_sel_genes(probes,find(ncount > 5));

   sel_arrays = 1:size(probes.pm,1);
   for i = 1:length(varargin)
      if(isstr(varargin{i}))
         switch(varargin{i})
            case 'steps',
               steps = varargin{i+1};
               i = i +1;
            case 'difference',
               difference = 1;
            case 'sel_arrays',
                i = i + 1;
                sel_arrays = varargin{i};
        end;
      end;
   end;

  
   narrays = size(probes.pm,1);
   ngenes = length(probes.name);

   % premalloc memory
   r = zeros(ngenes,narrays,steps);

   %interpolation location steps
   xsteps = 0:(1/(steps-1)):1;
  
   narray = size(probes.pm,1);

   signal = mt_real_signal(probes);
   
   narray = size(probes.pm,1);
   nprobe = size(probes.pm,2);
   if(~isfield(probes,'position'))
       warning('No position field, assuming probes are ordered according to 5-3')
   end;
    
   %filling data structure
   fprintf(1,' Creating data structure\n');
   for i = 1:length(probes.name)
      if ((mod(i,1000)==0)), fprintf(1,'.'); end;
      t = find(probes.ind == i);
      if(isfield(probes,'position'))
          w = sortrows([t probes.position(t)],2);
      else
          w = t(:);
      end;
      t = w(:,1);
      d = signal(:,t); 

      %determine if interpolation is necessary
      if(length(t) == steps)  
         r(i,:,:) = d;
      else
         w = 0:(1/(size(t,1) - 1)):1;
         r(i,:,:) = interp1(0:(1/(length(t)-1)):1,d',xsteps)';
      end;
   end;
   fprintf(1,'\n');

   res = squeeze(median(r,1));
   
   % if probes has labels for the arrays, use them to give the different 
   % array sets a different color
   d = squeeze(median(r,1));
   x = jet(length(sel_arrays));
   formats = {'-' ,':', '-.', '--'};
   if(exist('difference'))
       d = d - repmat(mean(d),narray,1);
   end;

   for a = 1:length(sel_arrays)
       i = sel_arrays(a);
       plot(d(i,:),formats{mod(i, 4) + 1}, 'Color',x(a,:));
       hold on;
   end;
   hold off;

   xlabel('Probe position order (5` <-- --> 3`)')
   ylabel('Log_2 expression intensity')
   title('Amplification effect')
   if(length(sel_arrays) < 30)
       legend(probes.array_names(sel_arrays), 'Location', 'EastOutside');
   end;

