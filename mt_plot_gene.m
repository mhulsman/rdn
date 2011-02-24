%MT_PLOT_GENE - Plots probeset and its summarized signal if available
%
%  MT_PLOT_GENE(PROBES,GENENR,VARARGIN)
%
% INPUT
%   PROBES			Probe structure
%   GENENR         Number of gene/probeset to visualize 
%   VARARGIN        'sort'     : Sort on summarized signal (or if not available, on probe with highest signal)
%                   'label'  :   If probeset contains annots cell array, with the various annotations per array
%                                in each cell (e.g. {{'M','W','M'},{'+','+','-'},[1930,1953,1980]}, then
%                                one can plot such labels with this option, e.g. mt_plot_gene(probes,100,'label', 1) 
%                                to plot first labelset.
%                                     
%
% DESCRIPTION
% Plots probeset log2 expression intensity for each probe, as well as the summarized
% signal if avialable. Probes sequences (if available) are shown in the legend. 
% 
% SEE ALSO
% MT_SUM_PLM, MT_PLOT_SEQ_EFFECTS, MT_PLOT_AMP_EFFECTS
%
% (c) Marc Hulsman, 2009
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function mt_plot_gene(probes,genenr,varargin)
   for i = 1:length(varargin)
      if(isstr(varargin{i}))
         switch(varargin{i})
            case 'sort',
               use_sort = 1;
            case 'label',
               i = i+1;
               use_label = varargin{i};
         end;
      end;
   end;

   narray = size(probes.pm,1);
   probe_idxs = find(probes.ind == genenr);
   cmap = jet(length(probe_idxs));

   signal = mt_real_signal(probes);
   signal = signal(:,probe_idxs);

   if(exist('use_sort'))
      if(isfield(probes,'array_factors'))
         m = probes.overall_factors(genenr) + probes.array_factors(:,genenr);
         [dummy,sidx] = sort(m);
         signal = signal(sidx,:);
      else
          [m,i] = max(max(signal));
          [signal,sidx] = sortrows(signal,i);
      end;
   else
      sidx = (1:narray)';
   end;
   
   style_pos = {'-','--','-.'};
   for i = 1:length(probe_idxs)
      style = style_pos{mod(i-1,length(style_pos)) + 1};
      h = plot(signal(:,i),style);
      set(h,'Color',cmap(i,:));
      hold on;
   end;

   if(isfield(probes,'array_factors'))
        m = probes.overall_factors(genenr) + probes.array_factors(:,genenr);
        m = m(sidx);
        h = plot(m,'k');
        set(h,'LineWidth',2);
   end;


   legend_ids = {};
   if(isfield(probes, 'sequence'))
    for i = 1:length(probe_idxs)
        legend_ids{i} = [num2str(i) ': ' probes.sequence{probe_idxs(i)}];
    end;
   else
    for i = 1:length(probe_idxs)
        legend_ids{i} = num2str(i);
    end;
   end;
   legend(legend_ids, 'Location', 'SouthOutside')
   grid on

   if(exist('use_label'))
      ylim = get(gca,'YLim');
      d = probes.annots(use_label);
      for i =1:length(sidx)
         text(i,20 + 0.5,d{sidx(i)});
         text(i,20 + 0.2,num2str(sidx(i)));
      end;
   end;

   xlabel('Array number')
   ylabel('Log_2 intensity')
   
   %set xtick labels according to sort order
   set(gca,'XTick',(1:narray));
   set(gca,'XTickLabel',num2str(sidx));

   %generating title
   if(isfield(probes, 'gene_names') & ~strcmp(probes.gene_names{genenr}, '---'))
       k = {[probes.gene_names{genenr},': ',probes.gene_description{genenr}]};

   else
    if(isfield(probes,'desc'))
        k = regexp(probes.desc{genenr},'DEF\=(.*?)\s/','tokens');
        if(length(k) == 0)
            k = regexp(probes.desc{genenr},'GEN\=(.*?)\s/','tokens');
        end;
        if(length(k) == 0)
            k = regexp(probes.desc{genenr},'UG_TITLE\=(.*?)\s/','tokens');
        end;
    else
        k = probes.name{genenr};
    end;
  end;

    if(length(k) == 1)
        title(['Gene plot of: ', k{1}])
    else
        title('Gene plot (gene unkown)')
    end;
   hold off;
   


