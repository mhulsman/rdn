%MT_SIG - Determines significant genes w.r.t. to some class labels
%
%  [BG,SSTAT] = MT_SIG_PERM(EXPRES,ANNOTNR,VARARGIN)
%
% INPUT
%   EXPRES			Expression/probes structure with summarized probeset signals
%   ANNOTNR         Annotation label nr in express structure. Probes requires annots field, with
%                   cell array containing label sets.
%                   For rank correlation and other tests: cell contains vector
%                   For 2-class tests (e.g. fold change): cell contains matrix (sample * 2) with 0/1 indicating membership of class
%   VARARGIN        'statfun'    :   Use different stat function (default: SAM test statistic (a regularized t-test))
%
% OUTPUT
%   BG		    Genes ordered on significance (top significant first)
%   SSTAT       Statistics for each probeset (ordered on BG)
%
% DESCRIPTION
%   Determines label significance for each gene. 
%
% (c) Marc Hulsman, 2010
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function [bg,sstat] = mt_sig(expres,annotnr,varargin)

%default stat function is sam statistic (a regularized t-test)
statfun = @mt_stat_sam;

%parse variables
for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'statfun',
            i = i + 1;
            statfun = varargin{i};
      end;
   end;
end;

%get labels
lbl = expres.annots{annotnr};
%get signal
signal = repmat(expres.overall_factors,size(expres.array_factors,1),1) + expres.array_factors;

%if 2 or more classes
if(size(lbl,2) >= 2)
    %remove samples without class
    not_remove = ~~sum(lbl,2);
    lbl = lbl(not_remove,:);

    %if one of classes is empty, return
    clabel = sum(lbl);
    if(clabel(1) == 0 | clabel(2) == 0)
        bg = [];
        pval = [];
        fdr = [];
        sstat = [];
        maxuniq = 0;
        return;
    end;
else
    %no classes, but continous.
    %remove nan samples
    not_remove = ~isnan(lbl);
    lbl = lbl(not_remove);
end;
signal = signal(not_remove,:);

[narray,ngene] = size(signal);

%calculate unpermuted stats
[sstat,xvar] = feval(statfun,signal,lbl,[]);

%sort over p-value and sstat
r = flipud(sortrows([(1:ngene)' abs(sstat) sstat],[2]));
bg = r(:,1);
sstat = r(:,3);

return;



