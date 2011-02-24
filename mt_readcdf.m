%MT_READCDF - Read an Affymetrix CDF file
%
%  CDF = MT_READCDF(FNAME,PROGRESS)
%
% INPUT
%   FNAME			File name 
%		PROGRESS	Show progress (1, default) or not (0)
%
% OUTPUT
%		CDF				CDF-structure
%
% DESCRIPTION
% Reads in an Affymetrix chip description (CDF) file. All fields are stored 
% in a self-explanatory structure. 
%
% NOTE
% Perfect Matches (PMs) are probes for which PBASE is the natural complement
% of TBASE in the CDF file. Mismatches (MMs) are probes for which this is not 
% true.
%
% SEE ALSO
% MT_READCEL, MT_CEL2PROBES

% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function cdf = mt_readcdf(fname,progress)

	if (nargin < 2), progress = 1; end;

	f = fopen(fname,'r');
	if (f <= 0), error ('cannot open file'); end;
	frewind(f);

    magic = fread(f, 1, 'int32');
    if(magic == 67)
        nblocks = 0;
        version = fread(f, 1, 'int32');
        if(version > 3)
            1/0
        end;
        cdf.cols = fread(f, 1, 'uint16');
        cdf.rows = fread(f, 1, 'uint16');
        cdf.nunits = fread(f, 1, 'int32');%excluding qc probesets
        cdf.qc_nunits = fread(f, 1, 'int32');
        custseq = readstr(f);
        probesetnames = cell(cdf.nunits, 1);
        for i = 1:cdf.nunits
            probesetnames{i} = char(remove_trail(fread(f, 64, 'uint8')))';
        end;
        qc_pos = fread(f, cdf.qc_nunits, 'int32');
        unit_pos = fread(f, cdf.nunits, 'int32');
        
        cdf.gene = cell(cdf.nunits, 1);

        for i = 1:cdf.nunits
            fseek(f, unit_pos(i), 'bof');
            cdf.gene{i}.name = probesetnames{i};

            utype = fread(f, 1, 'uint16');
            cdf.gene{i}.direction = fread(f, 1, 'uint8');
            cdf.gene{i}.natoms = fread(f, 1, 'int32');
            cdf.gene{i}.nblocks = fread(f, 1, 'int32');
            cdf.gene{i}.ncells = fread(f, 1, 'int32');
            cdf.gene{i}.unitno = fread(f, 1, 'int32');
            ncells_per_atom = fread(f, 1, 'uint8');
            cdf.gene{i}.block = cell(cdf.gene{i}.nblocks, 1);
            nblocks = nblocks + cdf.gene{i}.nblocks;

            for j = 1:cdf.gene{i}.nblocks
                cdf.gene{i}.block{j}.blockno = fread(f, 1, 'int32');
                cdf.gene{i}.block{j}.ncells = fread(f, 1, 'int32');
                cdf.gene{i}.block{j}.ncells_per_atom = fread(f, 1, 'uint8');
                cdf.gene{i}.block{j}.direction = fread(f, 1, 'uint8');
                first_atom_pos = fread(f, 1, 'int32');
                fread(f, 1, 'int32');
                cdf.gene{i}.block{j}.name = char(remove_trail(fread(f, 64, 'uint8')))';
                if(version > 1)
                    wobble = fread(f, 1, 'uint16');
                    allele_code = fread(f, 1, 'uint16');
                    if(version > 2)
                        channel = fread(f, 1, 'uint8');
                        reptype = fread(f, 1, 'uint8');
                    end;
                end;
                if(cdf.gene{i}.block{j}.ncells_per_atom > 1)
                    cdf.gene{i}.block{j}.pm = zeros(cdf.gene{i}.block{j}.ncells/2, 2,'uint16');
                    cdf.gene{i}.block{j}.mm = zeros(cdf.gene{i}.block{j}.ncells/2, 2,'uint16');
                else
                    cdf.gene{i}.block{j}.pm = zeros(cdf.gene{i}.block{j}.ncells, 2, 'uint16');
                end;

                for k = 1:cdf.gene{i}.block{j}.ncells
                    atom = fread(f, 1, 'int32');
                    x = fread(f, 1, 'uint16');
                    y = fread(f, 1, 'uint16');
                    atom_number2 = fread(f, 1, 'int32');
                    pbase = upper(char(fread(f, 1, 'uint8')));
                    tbase = upper(char(fread(f, 1, 'uint8')));
                    if(version >= 2)
                        probeseq_len = fread(f, 1, 'uint16');
                        physical_groupe = fread(f, 1, 'uint16');
                    end;
                    if (((pbase == 'A') & (tbase == 'T')) | ...
                        ((pbase == 'T') & (tbase == 'A')) | ...
                        ((pbase == 'G') & (tbase == 'C')) | ...
                        ((pbase == 'C') & (tbase == 'G')))
                        cdf.gene{i}.block{j}.pm(atom + 1,:) = [x+1 y+1];
                    else
                        pbase
                        tbase
                        cdf.gene{i}.block{j}.mm(atom + 1,:) = [x+1 y+1]
                    end;

                end;
            end;
        end;
        

    else
        frewind(f);
        % 1. Read version.
        if (search(f,'[CDF]') < 0), error ('cannot find [CDF] block'); end;
        line = fgets(f); cdf.version = sscanf(line,'Version=%s');
        ver = sscanf(cdf.version,'GC%f');
        
        % 2. Read chip info.
        
        if (search(f,'[Chip]') < 0), error ('cannot find [Chip] block'); end;
        line = fgets(f); cdf.name 								= sscanf(line,'Name=%s');
        line = fgets(f); cdf.rows 								= sscanf(line,'Rows=%d');
        line = fgets(f); cdf.cols 								= sscanf(line,'Cols=%d');	line = fgets(f); cdf.nunits 							= sscanf(line,'NumberOfUnits=%d');
        line = fgets(f); cdf.maxunit 							= sscanf(line,'MaxUnit=%d');
        line = fgets(f); cdf.qc_nunits 						= sscanf(line,'NumQCUnits=%d');
        line = fgets(f); cdf.reference 						= sscanf(line,'ChipReference=%d');

        % 3. Read QC units.

        for i = 1:cdf.qc_nunits
        if (search(f,'[QC') < 0), error ('cannot find next [QC] block'); end;
            % Skip for now.
        end;

        % 4. Read other units.

        if (progress), fprintf (1,'Reading units '); end;

        cdf.gene = cell(cdf.nunits,1);											% Allocate memory.
        
        maxatom = -1; nblocks = 0;

        for i = 1:cdf.nunits

            if (progress & (mod(i,1000)==0)), fprintf (1,'.'); end;

            if (search(f,'[Unit') < 0), error ('cannot find next [Unit] block'); end;
            line = fgets(f); cdf.gene{i}.name				= sscanf(line,'Name=%s');
            line = fgets(f); cdf.gene{i}.direction	= sscanf(line,'Direction=%d');
            line = fgets(f); cdf.gene{i}.natoms			= sscanf(line,'NumAtoms=%d');
            line = fgets(f); cdf.gene{i}.ncells			= sscanf(line,'NumCells=%d');
            line = fgets(f); cdf.gene{i}.unitno			= sscanf(line,'UnitNumber=%d');
            line = fgets(f); cdf.gene{i}.unittype		= sscanf(line,'UnitType=%d');
            line = fgets(f); cdf.gene{i}.nblocks		= sscanf(line,'NumberBlocks=%d');

            cdf.gene{i}.block = cell(cdf.gene{i}.nblocks,1);	% Allocate memory.

            nblocks = nblocks + cdf.gene{i}.nblocks;

            for j = 1:cdf.gene{i}.nblocks
                
            if (search(f,'[Unit') < 0), error ('cannot find next [Unit_Block] block'); end;
            line = fgets(f); cdf.gene{i}.block{j}.name 			= sscanf(line,'Name=%s'); 
                line = fgets(f); cdf.gene{i}.block{j}.blockno 	= sscanf(line,'BlockNumber=%d'); 
                line = fgets(f); cdf.gene{i}.block{j}.natoms 		= sscanf(line,'NumAtoms=%d'); 
                line = fgets(f); cdf.gene{i}.block{j}.ncells 		= sscanf(line,'NumCells=%d'); 
                line = fgets(f); cdf.gene{i}.block{j}.startpos 	= sscanf(line,'StartPosition=%d'); 
                line = fgets(f); cdf.gene{i}.block{j}.stoppos 	= sscanf(line,'StopPosition=%d'); 

                line = fgets(f);			% Skip line.

                cdf.gene{i}.block{j}.pm = zeros(cdf.gene{i}.block{j}.ncells/2,2,'uint16');
                cdf.gene{i}.block{j}.mm = zeros(cdf.gene{i}.block{j}.ncells/2,2,'uint16');

                atom_offset = 1;

                for k = 1:cdf.gene{i}.block{j}.ncells

                line = fgets(f); 
                    val  = sscanf(line,'Cell%d=%d %d %*c %*s %*s %d %d %c %c %c %d %d %d %d %d');
                    x = val(2); y = val(3); pos = val(4); 
                    pbase = val(6); tbase = val(7); atom = val(9);

                    % In version GC3.0. atom starts at 0; in GC2.0 it can start at 1.
                    % Only set the offset to 0 if the first atom starts at 1 (sigh)...
                    if ((ver < 3) & (k == 1) & (atom == 1)), atom_offset = 0; end;

                    atom = atom + atom_offset;

                    if (((pbase == 'A') & (tbase == 'T')) | ...
                        ((pbase == 'T') & (tbase == 'A')) | ...
                        ((pbase == 'G') & (tbase == 'C')) | ...
                        ((pbase == 'C') & (tbase == 'G')))
                        cdf.gene{i}.block{j}.pm(atom,:) = [x+1 y+1];
                    else
                        cdf.gene{i}.block{j}.mm(atom,:) = [x+1 y+1];
                    end;

                    maxatom = max(atom,maxatom);

                end;
            end;
        end;

        cdf.maxncells = maxatom;
        cdf.nblocks		= nblocks;
    end;
	fclose(f);

return

% SEARCH(F,STR)
% Look for the starting position of a string STR in file with handle F.
% Returns -1 if not found.

function pos = search (f, str)

	n = length(str);

	startpos = ftell(f);

	line = fgets(f); 
	while ((~strncmp(line,str,n)) & ~feof(f))
		line = fgets(f); 
	end;
	while ((~strncmp(line,str,n)) & (ftell(f) <= startpos))
		line = fgets(f); 
	end;

	if (strncmp(line,str,n)), pos = ftell(f); else, pos = -1; end;

return


function s = readstr(f)
    length = fread(f, 1, 'int32', 'ieee-be');
    s = char(fread(f, length, 'char'))';
return

function s = remove_trail(s)
    s = s(1:max(find(s ~= 0)));
return;
