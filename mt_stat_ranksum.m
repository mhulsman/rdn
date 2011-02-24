%MT_STAT_RANKSUM - Compute RANKSUM-statistics (Mann-Whitney U test)
%
%   S = MT_STAT_RANKSUM(signal,labels,s0)
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

      a1 = signal(~~labels(:,2),:); 
      a2 = signal(~~labels(:,1),:);
      [r,c] = size(a1);
      sval = zeros(c,1);

      for i = 1:c
        k1 = a1(:,i);
        k2 = a2(:,i);
        sval(i) = 1 - ranksum(k1,k2);
        if(mean(k2) < mean(k1))
            sval(i) = - sval(i);
        end;
      end;

