%MT_STAT_AUC - Compute AUC-statistics 
%
%   S = MT_STAT_AUC(signal,labels,s0)
%
% INPUT
%   signal  Array signal
%   labels  Class labels (2-class)
% 	s0      Extra variable
% OUTPUT
%   SVAL    RANKSUM-statistic for each probeset
%
% NOTE
%   slow. 
%
% SEE ALSO
% MT_STAT_SAM

% (c) Marc Hulsman, 2009
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function [sval,s0] = mt_stat_ranksum(signal,labels,s0)
	% Check input.
      s0 = [];
      signal = signal(~~sum(labels,2),:);
      labels = labels(~~sum(labels,2),:);
      
      w = sum(labels);
      poscount = w(1);
      negcount = w(2);
      labels = labels * [1;0];

      [dummy,i] = sort(signal);
      slabels = labels(i);

      tpsum = cumsum(slabels)/poscount;
      fpsum = cumsum(1-slabels)/negcount;

      dx = diff(fpsum);
      sval = sum(dx .* tpsum(1:(end-1),:))';

      low_filter = sval < 0.5;
      sval(low_filter) = 1 - sval(low_filter);

