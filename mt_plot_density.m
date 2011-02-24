%MT_PLOT_DENSITY - Plots signal density w.r.t intensity
%
%  MT_PLOT_CORR(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure
%   VARARGIN        %                                     
%
% DESCRIPTION
%    Plots density of signal w.r.t. expression intensity
%    Algorithm is a bit crude (histogram counting)
%
% SEE ALSO
% MT_COR_IMAGE, MT_PLOT_SEQ_EFFECTS, MT_PLOT_AMP_EFFECTS, MT_PLOT_GENE
%
% (c) Marc Hulsman, 2011
% Delft Bioinformatics Lab
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function res = mt_plot_density(probes)

signal = mt_real_signal(probes);

start = floor(min(min(signal)));
stop = ceil(max(max(signal)));

k = histc(signal',start:0.1:stop);

x = jet(size(signal,1));
formats = {'-' ,':', '-.', '--'};
for i = 1:size(signal,1)
    plot(k(:,i),formats{mod(i, 4) + 1}, 'Color',x(i,:));
    hold on;
end;
hold off;
xlabs = start:0.1:stop;
set(gca,'XTick',[1:10:length(xlabs)]);
set(gca,'XTickLabel',xlabs(1:10:length(xlabs)));
xlabel('Log2 expression')
ylabel('Probe count')
title('Signal density');

if(size(signal,1) < 30)
    legend(probes.array_names, 'Location', 'EastOutside');
end;

