%MT_PLOT_CORR - Plots correlation between arrays
%
%  MT_PLOT_CORR(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure
%   VARARGIN        %                                     
%
% DESCRIPTION
%    Plot correlation between probes based on (possibly normalized) signal. 
%
% SEE ALSO
% MT_COR_IMAGE, MT_PLOT_SEQ_EFFECTS, MT_PLOT_AMP_EFFECTS, MT_PLOT_GENE
%
% (c) Marc Hulsman, 2011
% Delft Bioinformatics Lab
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function res = mt_plot_corr(probes)

signal = mt_real_signal(probes);
res = corrcoef(signal');
imagesc(res);
colorbar
title('Correlation between arrays');
xlabel('Array nrs');
ylabel('Array nrs');
