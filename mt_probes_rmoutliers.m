%PROBES_RMOUTLIERS - Remove outliers
%
%   PROBES = PROBES_RMOUTLIERS (PROBES)
% 
% INPUT
%   PROBES		A probe structure containing PM and MM values (see CEL2PROBES)
%             and an OL array (see CEL2PROBES).
%
% OUTPUT
%		PROBES		Idem, containing normalised PM values
%
% DESCRIPTION
% Removes the outliers (indices stored in PROBES) between the probes.
%
% SEE ALSO
% CEL2PROBES, PROBES_NORM

% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function probes = probes_rmoutliers (probes)

	nprobes = length(probes.ind);

	if (isfield(probes,'ol') && ~isempty(probes.ol))
        probes.pm(:,probes.ol) = [];
		probes.ind(probes.ol)  = [];
        if(isfield(probes,'mm'))
	        probes.mm(:,probes.ol) = [];
        end;
		if (isfield(probes,'pm_pos'))
			probes.pm_pos(:,probes.ol) = [];
			probes.mm_pos(:,probes.ol) = [];
		end;
		if (isfield(probes,'position'))
			probes.position(probes.ol) = [];
		end;
		if (isfield(probes,'indmm'))
			probes.indmm(probes.ol) = [];
		end;
		if (isfield(probes,'sequence'))
			probes.sequence(probes.ol) = [];
		end;
		probes = rmfield(probes,'ol');
	else
		fprintf(1,'Outliers already removed...\n');
	end;
    
return

