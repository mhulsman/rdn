%MT_FIND_TAIL - determines probe distance to 3'end of probeset
%
%  PROBES = MT_FIND_TAIL(PROBES,VARARGIN)
%
% INPUT
%   PROBES			Probe structure including probe position
%
% OUTPUT
%   PROBES		    Input probe structure with result stored in lastpos
%
% DESCRIPTION
% Used to determine distance to poly-A tail. This version 
% simplifies the problem by assuming that this is close to the
% 3'end of the probeset
%
% SEE ALSO
% MT_CREATE_AMP_MODEL

% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands



function probes = mt_find_tail(probes);

nprobe = size(probes.pm,2);
t = sparse((1:nprobe)',probes.ind,(1:nprobe)');
last_probe = full(max(t));
probes.lastpos = zeros(size(probes.position));
for i = 1:length(probes.name)
   k = find(probes.ind == i);
   probes.lastpos(k) = double(probes.position(last_probe(i)) - probes.position(k));
end;
