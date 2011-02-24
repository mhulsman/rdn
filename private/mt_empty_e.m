%MT_EMPTY_E - Returns an (empty) expression structure from a probes structure
%
%  E = MT_EMPTY_E(PROBES,VARGIN)
%
% INPUT
%   PROBES		Probe structure used as basis
%
% OUTPUT
%   E		    Expression structure
%
% DESCRIPTION
% An expression structure is an probes structure with all probe level information
% removed
%
% SEE ALSO
% MT_SUM_PLM

% (c) Marc Hulsman, 2009
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function e = mt_empty_e(probes)

fn = fieldnames(probes);
remove = {'pm','mm','image_factors','qq_factors','seq_correction','amp_correction','seq_factors','seq_norm',...
    'seqbg','seqbg_factors','ind','lastpos','gene_sequence','bg_estimate','seqlength','pm_pos','mm_pos','sequence',...
    'position','nrows','ncols'};
remove = intersect(fn,remove);
e = rmfield(probes,remove);
