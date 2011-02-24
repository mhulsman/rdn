%MT_REAL_SIGNAL - Applies corrections and returns log2-transformed pm signal
%
%  YL = MT_REAL_SIGNAL(PROBES,VARGIN)
%
% INPUT
%   PROBES		Probe structure for which to calculate qq cor factors
%   VARARGIN    Optional string options, 'no_qq','no_image','no_sequence','no_amp'
%
% OUTPUT
%   YL		    Matrix (narray * nprobe) with corrected log2-transformed signal
%
% DESCRIPTION
% Calculates Quantile-Quantile correction factors (after application of ealier 
% corrections) and adds them to the probes structure as field qq_factors.   
%
% SEE ALSO
% MT_BG_EST, MT_COR_IMAGE, MT_COR_HYBAMP, MT_COR_QQ

% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands


function [yl] = mt_real_signal(probes,varargin)

nprobe = size(probes.pm,2);
narray = size(probes.pm,1);

for i = 1:length(varargin)
   if(isstr(varargin{i}))
      switch(varargin{i})
         case 'no_qq',
            no_qq = 1;
         case 'no_image',
            no_image = 1;
         case 'no_sequence',
            no_sequence = 1;
         case 'no_amp',
            no_amp = 1;
         case 'no_log',
            no_log = 1;
      end;
   end;
end;

if(exist('no_log'))
    yl = probes.pm;
else
    yl = log2(probes.pm);
end;

if(isfield(probes,'qq_factors') && ~exist('no_qq'))
    fprintf(1,' Correcting for distribution effects\n');
    yl = yl - probes.qq_factors;
end;

if(isfield(probes,'image_factors') && ~exist('no_image'))
    fprintf(1,' Correcting for image effects\n');
    yl = yl - probes.image_factors;
end;

if(isfield(probes,'ampseq_correction') && ~exist('no_sequence') && ~exist('no_amp'))
    fprintf(1,' Performing sequence-based correction\n');
    yl = yl - probes.ampseq_correction;
end;

if(isfield(probes,'seq_correction') && ~exist('no_sequence'))
    fprintf(1,' Performing sequence-based correction\n');
    yl = yl - probes.seq_correction;
end;
if(isfield(probes,'amp_correction') && ~exist('no_amp'))
    fprintf(1,' Performing amplification-based correction\n');
    yl = yl - probes.amp_correction;
end;

