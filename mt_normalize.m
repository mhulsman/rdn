%MT_NORMALIZE - performs all normalization steps
%
%  PROBES = MT_NORMALIZE(PROBES);
%
% INPUT
%   PROBES			Probe structure for which to calculate qq cor factors
% OUTPUT
%   PROBES		    Normalized probes structure
%
% SEE ALSO
% MT_COR_QQ, MT_COR_IMAGE, MT_COR_HYBAMP, MT_REAL_SIGNAL

% (c) Marc Hulsman, 2010
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function probes = mt_normalize(probes, varargin)


for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'no_summarization',
            no_summarization = 1;
      end;
   end;
end;

probes = mt_bg_est(probes,'reg_max',10,'m_estimation');
if isfield(probes, 'gene_sequence')
    probes = mt_cor_hybamp(probes,'m_estimation','use_amplification');
else
    probes = mt_cor_hybamp(probes,'m_estimation');
end;
probes = mt_cor_image(probes,9,'use_optical');
probes = mt_cor_qq(probes);

if ~exist('no_summarization')
    probes = mt_sum_plm(probes,'m_estimation','keep_probe_info');
end;    
