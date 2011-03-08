%MT_EXPR_SIGNAL - Returns summarized signal
%
%  YL = MT_EXPR_SIGNAL(PROBES)
%
% INPUT
%   PROBES		Probe structure for which to calculate expression signal
%
% OUTPUT
%   YL		    Matrix (narray * ngene) with summarized signal
%
% (c) Marc Hulsman, 2011
% Delft Bioinformatics Lab
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands


function [yl] = mt_expr_signal(probes)


yl = probes.array_factors + repmat(probes.overall_factors, size(probes.array_factors, 1), 1);

