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
   order = 1:length(probes.array_names);

   for i = 1:length(varargin)
      if(isstr(varargin{i}))
         switch(varargin{i})
            case 'sort',
               use_sort = 1;
            case 'label',
               i = i+1;
               use_label = varargin{i};
            case 'order',
                i = i+1;
                lbls = probes.annots{varargin{i}};
                if(size(lbls,2) > 1)
                    lbls = lbls * (1:size(lbls,2))';
                end;
                [dummy, order] = sort(lbls);
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
         [dummy,order] = sort(m);
      end;
   end;

   
   
   style_pos = {'-','--','-.'};
   for i = 1:length(probe_idxs)
      style = style_pos{mod(i-1,length(style_pos)) + 1};
      h = plot(signal(order,i),style);
      set(h,'Color',cmap(i,:));
      hold on;
   end;

   if(isfield(probes,'array_factors'))
        m = probes.overall_factors(genenr) + probes.array_factors(order,genenr);
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


   xlabel('Array number')
   ylabel('Log_2 intensity')
   
   %set xtick labels according to sort order
   jump = ceil(narray / 40);
   set(gca,'XTick',(1:jump:narray));
   if(exist('use_label'))
      lbl = probes.annots{use_label};
      if(size(lbl,2) > 1)
        lbl = lbl * (1:size(lbl,2))';
      end;
      xlabs = num2str(lbl(order));
   else
      xlabs = num2str((1:narray)');
   end;
   xlabs = xlabs(1:jump:narray,:);
   set(gca,'XTickLabel',xlabs);
   

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
            if(length(k) == 0)
                k = {sprintf('%s: %s',probes.name{genenr},probes.desc{genenr})};
            end;
        else
            k = {sprintf('%s: %s',probes.name{genenr},probes.desc{genenr})};
        end;
   end;
   if(length(k) == 0)
       k= probes.name(genenr);                
   end;

   if(length(k) == 1)
       k{1} = strrep(k{1},sprintf('\t'),' ');
       title(['Gene plot of: ', k{1}])
   else
       title('Gene plot (gene unkown)')
   end;
   hold off;
   


