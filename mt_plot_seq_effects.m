%MT_PLOT_SEQ_EFFECTS - Plots sequence effects for all arrays
%
%  MT_PLOT_SEQ_EFFECTS(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure for which to calculate cor factors
%   VARARGIN        'chars'     : nucelotides to count ('A','C', 'G',  and/or 'T').
%                                 default [2,3]  (i.e. 'C' and 'G')
%                   'difference'  : Show difference values w.r.t mean instead of absolute values   
%                   'seq_range' : Use only part of probes sequences (default: 1:25)
%                   'limit'     : If for certain nucelotide count, probe count < limit, then
%                                 do not show that nucleotide count. Default: 100
%                                     
%
% DESCRIPTION
%   Displays single-nucleotide sequence effects for probes with a 
%   specific amount of certain nucleotides.
%
% SEE ALSO
% MT_COR_HYBAMP, MT_PLOT_ARRAY, MT_PLOT_AMP_EFFECTS, MT_PLOT_GENE
%
% (c) Marc Hulsman, 2009
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function mt_plot_seqeffects(probes,varargin)

chars = [2,3];
seq_range = 1:25;
limit = 100;
for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'chars',
            chars = varargin{i+1};
            i = i + 1;
         case 'difference',
            difference = 1;
         case 'seq_range',
            i = i + 1;
            seq_range = varargin{i};
         case 'limit',
            i = i + 1;
            limit = varargin{i};
      end;
  end;
end;

dna_letters = ['A','C','G','T'];
q = upper(char(probes.sequence));

q = q(:,seq_range);

nprobe = size(probes.pm,2);
nucleo_count = zeros(nprobe,4);

for i = 1:length(dna_letters)
   t =  (q == dna_letters(i));
   nucleo_count(:,i) = sum(t,2);
end;

if(length(chars) > 1)
   gc_content = sum(nucleo_count(:,chars),2);
else
   gc_content = nucleo_count(:,chars);
end;

gc_u = unique(gc_content)';

gc_nr = histc(gc_content,gc_u)';
[gc_u; gc_nr]

remove = gc_nr < limit;
gc_u = gc_u(~remove);
gc_nr = gc_nr(~remove);

narray = size(probes.pm,1);

signal = mt_real_signal(probes);

data = zeros(narray,length(gc_u));

for i = 1:length(gc_u)
   data(:,i) = median(signal(:,find(gc_content == gc_u(i))),2);
end;

x = jet(size(data,1));
formats = {'-', ':', '-.', '--'};

if(exist('difference'))
    data = data - repmat(mean(data),narray,1);
end;
for i = 1:size(data,1)
    plot(data(i,:),formats{mod(i, 4) + 1},'Color',x(i,:));
    hold on;
end;
hold off;

set(gca,'XTick',1:length(gc_u))
set(gca,'XTickLabel',gc_u)
xlabel('Number of specified nucleotides in probe')
ylabel('Log_2 median expression intensity')
title(strcat('Sequence effect for nucleotide(s): ', dna_letters(chars)));
legend(probes.array_names, 'Location','EastOutside')

