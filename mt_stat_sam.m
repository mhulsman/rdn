%MT_STAT_RANKSUM - Compute RANKSUM-statistics
%
%   S = MT_STAT_SAM(signal,labels,s0)
%
% INPUT
%   signal  Array signal
%   labels  Class labels (2-class)
% 	s0      Extra variable
%
% OUTPUT
%   SVAL    RANKSUM-statistic for each probeset
%
% DESCRIPTIOn
% S-statistics are the statistics used in SAM. Basically, they're two-tailed
% pooled-variance T-statistics, but with an estimated "fudge factor" S0
% in the denominator.
%
% SEE ALSO
% MT_STAT_RANKSUM
%
% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function [sval,s0] = mt_stat_sam(signal,labels,s0)

	% Check input.

      a1 = signal(~~labels(:,2),:); 
      a2 = signal(~~labels(:,1),:); 
      n = [size(a1,1),size(a2,1)];
      m1 = mean(a1,1); m2 = mean(a2,1);
      r = m2 - m1;
      s = sqrt ((((1/n(1))+(1/n(2))) / (n(1) + n(2) - 2)) * ...
                  (sum((a1-ones(n(1),1)*m1).^2) + ...
                  sum((a2-ones(n(2),1)*m2).^2)));
   
      if (isempty(s0))
         s0 = find_s0(r,s);
      end;
      
      sval = - (r' ./ (s'+s0));
      sval(isnan(sval)) = 0;

return

function s0 = find_s0(r,s)

	for alpha = 1:101, q(alpha) = prctile(s,alpha-1); end;
	
	alpha = 0:0.05:1;

	for a = 1:length(alpha)
		d = r./(s+q(floor(alpha(a)*100+1))+realmin);
		for j = 1:100
			ind1 = find(s>=q(j)); ind2 = find(s(ind1)<q(j+1)); ind = ind1(ind2);
			if (~isempty(ind))
				md   = median(d(ind)); 
				v(j) = median(abs(d(ind)-md)) / .64;
			else
				v(j) = 0;
			end;
		end;
		cv(a) = (std(v)+realmin)./(mean(v)+realmin);
	end;

	[dummy,ind] = min(cv); 
	s0 = prctile(s,floor(alpha(ind)*100));
   
return
