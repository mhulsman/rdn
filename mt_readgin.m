%MT_READGIN - Read an Affymetrix GIN file
%
%  GIN = MT_READGIN(FNAME,PROGRESS)
%
% INPUT
%   FNAME			File name 
%		PROGRESS	Show progress (1, default) or not (0)
%
% OUTPUT
%		GIN				Structure containing gene descriptions
%
% DESCRIPTION
% Reads in an Affymetrix gene description (GIN) file. All fields are stored 
% in a self-explanatory structure. 
%
% SEE ALSO
% MT_READCDF, MT_READCEL, MT_READPROBE_ANNOT, MT_CEL2PROBES

% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function gin = mt_readgin(fname,progress)

	if (nargin < 2), progress = 1; end;

	f = fopen(fname,'r');
	if (f <= 0), error ('cannot open file'); end;
	frewind(f);

	% Read version and skip next couple of lines.

	line = fgets(f); gin.version = sscanf(line,'Version=%s');
	while (~strncmp(line,'Index',5))
         line = fgets(f);
   end; 

	done = 0; 
	while (~feof(f))

		line = fgets(f);
		line = strrep(line,sprintf('\t\t'),sprintf('\tX\t'));		% Prevent empty fields.

		[n,cnt,err,ni] = sscanf(line, '%d ', 1); line = line(ni:end);

		[dummy,cnt,err,ni] = sscanf(line, '%s ', 1); line = line(ni:end);
		[dummy,cnt,err,ni] = sscanf(line, '%s ', 1); line = line(ni:end);

		[gin.name{n},cnt,err,ni] = sscanf(line, '%s ', 1); line = line(ni:end);
		line = strrep(line,sprintf('\tX\t'),sprintf('\t\t'));	% Replace empty fields.
		line(find(line==10)) = []; line(find(line==13)) = []; % Remove newlines
		gin.desc{n} = line;

	end;

	fclose(f);

return

