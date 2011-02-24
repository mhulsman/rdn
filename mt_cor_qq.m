%MT_COR_QQ - Calculates Quantile-Quantile correction factors
%
%  PROBES = MT_COR_QQ(PROBES)
%
% INPUT
%   PROBES			Probe structure for which to calculate qq cor factors
%
% OUTPUT
%   PROBES		    Probes structure with qq cor. factors added
%
% DESCRIPTION
% Calculates Quantile-Quantile correction factors (after application of ealier 
% corrections) and adds them to the probes structure as field qq_factors.   
%
% SEE ALSO
% MT_BG_EST, MT_COR_IMAGE, MT_COR_HYBAMP, MT_REAL_SIGNAL
% Revision history
%   Marc Hulsman - 8/2008 - Adapt for use with mt_real_signal

% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function probes = mt_cor_qq(probes,varargin);
    
    if(~isfield(probes,'qq_factors'))
      probes.qq_factors = zeros(size(probes.pm));
    end;

    yl = mt_real_signal(probes);
    
    [narrays,nprobes] = size(yl);

    [ssort,ind] = sort(yl');
    newsig = ones(narrays,1)*mean(ssort');
    new_yl = zeros(size(yl));
    for c = 1:narrays
        new_yl(c,ind(:,c)) = newsig(c,:);
    end;
    
    probes.qq_factors = probes.qq_factors + (yl - new_yl);


