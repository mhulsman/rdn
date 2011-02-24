%MT_SIG_PERM - Determines significant genes w.r.t. to some class labels using permutation
%
%  [BG,PVAL,FDR,SSTAT] = MT_SIG_PERM(EXPRES,NPERMS,ANNOTNR,VARARGIN)
%
% INPUT
%   EXPRES			Expression/probes structure with summarized probeset signals
%   NPERMS          Number of permutations
%   ANNOTNR         Annotation label nr in express structure. Probes requires annots field, with
%                   cell array containing label sets.
%                   For rank correlation and other tests: cell contains vector
%                   For 2-class tests (e.g. fold change): cell contains matrix (sample * 2) with 0/1 indicating membership of class
%   VARARGIN        'statfun'    :   Use different stat function (default: SAM test statistic (a regularized t-test))
%
% OUTPUT
%   BG		    Genes ordered on significance (top significant first)
%   PVAL        P-value determined by permutation (ordered on BG)
%   FDR         False discovery rate (FDR) (ordered on BG)
%   SSTAT       Statistics for each probeset (ordered on BG)
%
% DESCRIPTION
%   Determines label significance for each gene using label permutation. Can scale to high number of permutations
%   by only using required number of permutations per gene to determine significance (i.e. 25 permutation tests with higher score).
%
% (c) Marc Hulsman, 2010
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function [bg,pval,fdr,sstat,maxuniq] = mt_sig_perm(expres,nperms,annotnr,varargin)

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


%determine maximal number of unique permutations

%3-classes (permute over all samples, but calculate score over class 1 and 2)
if(size(lbl,2) == 3)
        clabel = sum(lbl);
        maxuniq = factorial(size(lbl,1)) / (factorial(clabel(1)) * factorial(clabel(2)) * factorial(clabel(3)));
else 
    %2 classes
    if(size(lbl,2) == 2)
        clabel = sum(lbl);
        maxuniq = factorial(size(lbl,1)) / (factorial(clabel(1)) * factorial(clabel(2)));
    else
        %continous values
        maxuniq = factorial(size(lbl,1));
    end;
end;
disp(sprintf('Maximum number of unique permutations: %d',maxuniq));

%determine if number of requested permutations is larger than maximum number of unique permutations
if(maxuniq < nperms)
    warning('Number of permuatations %d larger than maximum number of unique permutations %d. Method will still work.',nperms,maxuniq);
    %more than 10 times is not necessary
    if(nperms > 10 * maxuniq)
        warning('Lowering max permutation count to %d',ceil(10 * maxuniq));
        nperms = ceil(10 * maxuniq);
    end;
end;


[narray,ngene] = size(signal);

%calculate unpermuted stats
[sstat,xvar] = feval(statfun,signal,lbl,[]);

%result vector with number of permutations better scoring than original
permres = zeros(length(sstat),1);

%we use an optimized algorithm, leaving out samples for which we already have enough permutations
%here, we use 25 permutations better than original as a threshold for that. 

%create signal indices for which we still have to perform permutations
filter = permres < 25;
fsignal = signal(:,filter);
fsstat = sstat(filter);

%run permutations
mlength = fprintf(1,'Permutating: %d/%d with %d probesets left...',0,nperms,sum(filter));
fpermres = zeros(size(permres));
for i = 1:nperms
    %every 1000 permutations, refilter the dataset
    if(mod(i,1000) == 0)
        %calculate new permres
        permres(filter) = permres(filter) + fpermres;
        %filter on threshold
        nfilter = permres < 25;
        %extrapolate previously filtered results
        extrap = filter & ~nfilter;
        permres(extrap) = (permres(extrap) / (i-1)) * nperms;
        %update filter
        filter = nfilter;
        nprobesets = sum(filter);
        fprintf(1,repmat('\b',1,mlength));
        mlength = fprintf(1,'Permutating: %d/%d with %d probesets left...',i,nperms,nprobesets);

        fpermres= zeros(nprobesets,1);
        %if no probesets left, we are done.
        if(nprobesets == 0)
            break;
        end;
        fsignal = signal(:,filter);
        fsstat = sstat(filter);
    end;
    pstat = feval(statfun,fsignal,lbl(randperm(narray),:),xvar);
    fpermres = fpermres + (abs(pstat) > abs(fsstat));
end;
%final permresults
permres(filter) = permres(filter) + fpermres;

%calculate p-values from permres
if(maxuniq < nperms)
    %correct for number of unique permutations
    permres = (permres / nperms) + (1 / (maxuniq + 1));
else
    permres = (permres + 1) / (nperms + 1);
end;
fprintf(1,'\n');

%sort on p-value
[pval,bg] = sort(permres);
sstat = sstat(bg);
%sort over p-value and sstat
r = flipud(sortrows([bg 1-pval abs(sstat) sstat],[2,3]));
bg = r(:,1);
pval = 1 - r(:,2);
sstat = r(:,4);
%calculate fdr
fdr = pval ./ ((1:ngene)' ./ ngene);
fdr = flipud(mt_cummin(flipud(fdr)));

return;



