%MT_READCEL - Read one or more Affymetrix CEL files
%
%  CEL = MT_READCEL(FNAME,PROGRESS)
%
% INPUT
%   FNAME			File name or cell array of file names
%		PROGRESS	Show progress (1, default) or not (0)
%
% OUTPUT
%		CEL				Cell array of CEL-structures
%
% DESCRIPTION
% Reads in a number of Affymetrix CEL files. All fields are stored in a
% self-explanatory structure.
%
% SEE ALSO
% MT_READCDF, MT_CEL2PROBES
%
% RVISION HISTORY
% M1, Marc, 12/2008 - Speedup
% M2, Michel Bellis, 4/2010 - Chip format detection
%
% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function cel = mt_readcel(fname,progress)

	if (nargin < 2), progress = 1; end;

	if (~iscell(fname)), fname = {fname}; end;

	for fn = 1:length(fname)

      	f = fopen(fname{fn},'r');
  	    if (f <= 0), error ('cannot open file'); end;
      	frewind(f);

      	if (progress), fprintf(1,'Array %d:\n',fn); end;
      	if (progress), fprintf(1,'   headers\n'); end;

      	% 1. Read version.
        magic = fread(f,1,'uint8');
        cel{fn}.version = fread(f, 1, 'uint8');
        if ((magic == 59) && (cel{fn}.version == 1))
            ndatagroup = fread(f, 1, 'int32', 'ieee-be');
            datagroup_pos = fread(f, 1, 'uint32', 'ieee-be');
            header = readheader(f);
            cel{fn}.header = header;

            %go to first data group (should already be there after reading
            %header, but just in case
            fseek(f, datagroup_pos, 'bof');
            group = readgroup(f);

            %somehow cel_cols and cel_rows do not give the correct answer...
            %try to find using Grid* params
            if(isfield(header,'algorithm_param_GridURX'))
                chipcols= double(round(header.algorithm_param_GridURX(1) - header.algorithm_param_GridULX(1)));
                chiprows = double(round(header.algorithm_param_GridLLY(1) - header.algorithm_param_GridULY(1)));
            else
                chipcols = 1;
                chiprows = 1;
            end;

            %M2
            %if chipcols*chiprows~=nb of probes => sqr(nb of probes)xsqr(nb of probes) is suggested as the right chip size
            %all chips must have the same size that the first chip (either calculated (chipcols, chiprows) or changed by user)
            if chiprows*chipcols~=length(group.datasets{1}.cols.data)
                if fn==1
                    currrows=round(sqrt(length(group.datasets{1}.cols.data)));
                    currcols=currrows;
                    if currrows*currcols==length(group.datasets{1}.cols.data)
                        warning('Dimensions (rows=%u * cols=%u) found by algoritm.', currrows, currcols);
                        %choice=questdlg(sprintf('do you accept these calculated dimensions rows=cols=%u instead of (rows=%u * cols=%u ~= %u found by algorithm)',currcols,chiprows,chipcols,length(group.datasets{1}.cols.data)),'','yes','no','cancel','yes');
                        %if isequal(choice,'cancel')
                        %    h=errodlg('process canceled');
                        %    waitfor(h)
                        %%    error('process canceled')
                        %elseif isequal(choice,'no')
                        %    [rows,cols]=enterdim(length(group.datasets{1}.cols.data));
                        %    currcols=cols;
                        %    currrows=rows;                           
                        %elseif isequal(choice,'yes')
                        cols=currcols;
                        rows=currrows;
                        %end
                    else
                        h=warndlg(sprintf('calculated row*cols (%u*%u=%u) ~= %u',cols,cols,cols*cols,itemb));
                        waitfor(h)
                        [currrows,currcols]=enterdim(length(group.datasets{1}.cols.data))
                        cols=currcols;
                        rows=currrows;
                    end
                elseif cols==0
                    %previous chip had different size
                    h=errordlg('dimensions are changing');
                    waitfor(h)
                    error('process canceled')
                else
                    currcols=cols;
                    currrows=rows;                           
                end                
            else
                if cols==0
                    currcols=chipcols;
                    currrows=chiprows;                           
                else
                    %previous chip had different size
                    h=errordlg('dimensions are changing');
                    waitfor(h)
                    error('process canceled')
                end
            end
                                
            cel{fn}.cols=currcols;
            cel{fn}.rows=currrows;           

        else
            fseek(f,0,'bof');
            magic = fread(f,1,'int32');
            cel{fn}.version = fread(f,1,'int32');
            if ((magic == 64) && (cel{fn}.version >= 4))
                % Read sizes.
                cel{fn}.cols        = fread(f,1,'int32');
                cel{fn}.rows        = fread(f,1,'int32');
                cel{fn}.cells       = fread(f,1,'int32');
                cel{fn}.header_len  = fread(f,1,'int32');
            else
                % Seek CEL block.
                fseek(f,0,'bof'); if (search(f,'[CEL]') < 0), error ('cannot find [CEL] block'); end;
                line = fgets(f); cel{fn}.version = sscanf(line,'Version=%d');
                % Seek header.
                if (search(f,'[HEADER]') < 0), error ('cannot find [HEADER] block'); end;
                
                % 2. Read header

                line = fgets(f); cel{fn}.cols 				    = sscanf(line,'Cols=%d');
                line = fgets(f); cel{fn}.rows 					= sscanf(line,'Rows=%d');
                line = fgets(f); cel{fn}.totalx 				= sscanf(line,'TotalX=%d');
                line = fgets(f); cel{fn}.totaly 				= sscanf(line,'TotalY=%d');
                line = fgets(f); cel{fn}.offsetx 				= sscanf(line,'OffsetX=%d');
                line = fgets(f); cel{fn}.offsety 				= sscanf(line,'OffsetY=%d');
                line = fgets(f); cel{fn}.gridcornerul 			= sscanf(line,'GridCornerUL=%d %d');
                line = fgets(f); cel{fn}.gridcornerur			= sscanf(line,'GridCornerUR=%d %d');
                line = fgets(f); cel{fn}.gridcornerlr 			= sscanf(line,'GridCornerLR=%d %d');
                line = fgets(f); cel{fn}.gridcornerll 			= sscanf(line,'GridCornerLL=%d %d');
                line = fgets(f); cel{fn}.axisinvertx 			= sscanf(line,'Axis-invertX=%d');
                line = fgets(f); cel{fn}.axisinverty 			= sscanf(line,'AxisInvertY=%d');
                line = fgets(f); cel{fn}.swapxy 				= sscanf(line,'swapXY=%d');
                line = fgets(f); cel{fn}.datheader 				= sscanf(line,'DatHeader=%s');
                n = regexp(line(11:end),'\s(?<name>[^\s]+):CLS','names');
                cel{fn}.name = n.name;
                cel{fn}.filename = fname{fn};
                line = fgets(f); cel{fn}.algorithm 				= sscanf(line,'Algorithm=%s');
                line = fgets(f); cel{fn}.algorithmparameters 	= sscanf(line,'AlgorithmParameters=%s');
            end;
        end;
        cel{fn}.name = fname{fn};
        cel{fn}.filename = fname{fn};

      	% Allocate memory.

      	cel{fn}.nprobes 	= cel{fn}.cols * cel{fn}.rows;

        switch(cel{fn}.version)
            case 1
                cel{fn}.noutliers = group.datasets{4}.nrows;
                cel{fn}.nmasks = group.datasets{5}.nrows;
                cel{fn}.mean = single(reshape(group.datasets{1}.cols.data, cel{fn}.rows, cel{fn}.cols));
                %cel{fn}.std = reshape(group.datasets{2}.cols.data, cel{fn}.rows, cel{fn}.cols));
                %cel{fn}.npixels = reshape(group.datasets{3}.cols.data, cel{fn}.rows, cel{fn}.cols);
          	    %cel{fn}.mask = spalloc(cel{fn}.rows,cel{fn}.cols,cel{fn}.nmasks);
                %pos = group.datasets{4}.cols;
                %for i = 1:cel{fn}.nmasks
                %    cel{fn}.mask(pos(1).data(i) + 1, pos(2).data(i) + 1) = 1;
                %end;
          	    cel{fn}.outlier = spalloc(cel{fn}.rows,cel{fn}.cols,cel{fn}.noutliers);
                pos = group.datasets{4}.cols;
              
                for i = 1:cel{fn}.noutliers
                    cel{fn}.outlier(pos(1).data(i) + 1, pos(2).data(i) + 1) = 1;
                end;

            case 4
                % Skip algorithm and parameter settings blocks.
                fseek(f,24+cel{fn}.header_len,'bof');  alg_len = fread(f,1,'int32');
                fseek(f,alg_len,'cof');        par_len = fread(f,1,'int32');
                fseek(f,par_len,'cof');        cell_margin = fread(f,1,'int32');

                % Read data sizes.
                cel{fn}.noutliers   = fread(f,1,'uint32');
                cel{fn}.nmasks      = fread(f,1,'uint32');
                cel{fn}.nsubgrids   = fread(f,1,'int32');
            
                if (progress), fprintf(1,'   %d probes ',cel{fn}.nprobes); end;

                % 3. Read intensities. Use the 'skip' option.

                pos = ftell(f);
                cel{fn}.mean = single(reshape(fread(f,cel{fn}.nprobes,'float32',6),cel{fn}.rows,cel{fn}.cols));
                %fseek(f,pos+4,'bof');
                %cel{fn}.std  = reshape(fread(f,cel{fn}.nprobes,'float32',6),cel{fn}.rows,cel{fn}.cols);
                %fseek(f,pos+8,'bof');
                %cel{fn}.npixels = reshape(fread(f,cel{fn}.nprobes,'int16',6),cel{fn}.rows,cel{fn}.cols);
                fseek(f,pos+10*cel{fn}.nprobes,'bof');

                if (progress), fprintf(1,'\n  masks, outliers, modified\n'); end;

                % 4. Read masked probes. Allocate sparse matrix.
                

                if (cel{fn}.nmasks > 0)
                    for i = 1:cel{fn}.nmasks
                        pos = fread(f,2,'int16');
                        cel{fn}.mask(pos(1)+1,pos(2)+1) = 1; 
                    end
                end;

                % 5. Read outliers. Allocate sparse matrix.
                
                cel{fn}.outlier = spalloc(cel{fn}.rows,cel{fn}.cols,cel{fn}.noutliers);	

                if (cel{fn}.noutliers > 0)
                    for i = 1:cel{fn}.noutliers
                        pos = fread(f,2,'int16');
                        cel{fn}.outlier(pos(1)+1,pos(2)+1) = 1; 
                    end
                end;        
            
            case 3
                cel{fn}.mean    	= zeros(cel{fn}.rows,cel{fn}.cols, 'single');
                %cel{fn}.std     	= zeros(cel{fn}.rows,cel{fn}.cols);
                %cel{fn}.npixels 	= zeros(cel{fn}.rows,cel{fn}.cols);

                % 3. Read intensities.

                if (search(f,'[INTENSITY]') < 0), error ('cannot find [INTENSITY] block'); end;
                line = fgets(f); 
                if (sscanf(line,'NumberCells=%d') ~= cel{fn}.nprobes)
                    error ('number of probes does not match');
                end;
                line = fgets(f);	% Skip legend line.
                
                if (progress), fprintf(1,'   %d probes ',cel{fn}.nprobes); end;
        
                %M1
                C = textscan(f,'%d %d %f %f %d');
                p1 = C{1}+1; p2 = C{2}+1; mev = C{3}; sev = C{4}; npix = C{5};
                pos = sub2ind(size(cel{fn}.mean),p1,p2);

                cel{fn}.mean(pos) 		= mev;		% Store values.
                %cel{fn}.std(pos)  		= sev;
                %cel{fn}.npixels(pos) 	= npix;
            
                if (progress), fprintf(1,'\n  masks, outliers, modified\n'); end;
                        
                % 4. Read masks.
            
                if (search(f,'[MASKS]') < 0), error ('cannot find [MASKS] block'); end;
                line = fgets(f); 
                cel{fn}.nmasks = sscanf(line,'NumberCells=%d');
            
                % Allocate sparse matrix.
                %cel{fn}.mask = spalloc(cel{fn}.rows,cel{fn}.cols,cel{fn}.nmasks);	
        
                if (cel{fn}.nmasks > 0)
                    line = fgets(f);	% Skip legend line.
                    for i = 1:cel{fn}.nmasks
                        line = fgets(f); 
                        val	 = sscanf(line,' %d %d');									% Read position.
                        %cel{fn}.mask(val(1)+1,val(2)+1) = 1;					% Mark as masked.
                    end;
                end;
        
                % 5. Read outliers.

                if (search(f,'[OUTLIERS]') < 0), error ('cannot find [OUTLIERS] block'); end;
                line = fgets(f); 
                cel{fn}.noutliers = sscanf(line,'NumberCells=%d');
        
                % Allocate sparse matrix.
                cel{fn}.outlier = spalloc(cel{fn}.rows,cel{fn}.cols,cel{fn}.noutliers);	

                if (cel{fn}.noutliers > 0)
                    line = fgets(f);	% Skip legend line.
                    for i = 1:cel{fn}.noutliers
                        line = fgets(f);
                        val  = sscanf(line,' %d %d');									% Read position.
                        cel{fn}.outlier(val(1)+1,val(2)+1) = 1;				% Mark as outlier.
                    end;
                end;
                
                % 6. Read modified.

                if (search(f,'[MODIFIED]') < 0), error ('cannot find [MODIFIED] block'); end;
                line = fgets(f); 
                cel{fn}.nmodified = sscanf(line,'NumberCells=%d');

                % Allocate sparse matrix.
                %cel{fn}.modified = spalloc(cel{fn}.rows,cel{fn}.cols,cel{fn}.nmodified);	

                if (cel{fn}.nmodified > 0)
                    line = fgets(f);	% Skip legend line.
                    for i = 1:cel{fn}.nmodified
                        line = fgets(f);
                        val  = sscanf(line,' %d %d %f');							% Read position.
                        %cel{fn}.modified(val(1)+1,val(2)+1) = val(3);	% Store original mean.
                    end;
                end;

                fclose(f);

            end;
    end;
    
    % If just one filename is specified, don't return a cell array of CELs.

	if (length(fname)==1)
		cel = cel{1};
	end;

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


function header = readheader(f)
    header.type = readstr(f);
    header.guid = readstr(f);
    header.creationdate = readwstr(f);
    header.locale = readwstr(f);

    %read params
    header = readparams(f,  header);
    if(isfield(header,'error_occurred'))
        header.parents = {};
        return;
    end;
    nparentheader = fread(f, 1, 'int32', 'ieee-be');
    parent = {};
    for i = 1:nparentheader
        parent{end + 3} = readheader(f);
    end;
    header.parents = parent;
return

function group = readgroup(f)
    group.nextgroup = fread(f, 1, 'uint32', 'ieee-be');
    dataset_pos = fread(f, 1, 'uint32', 'ieee-be');
    numdataset = fread(f, 1, 'int32', 'ieee-be');
    group.group_name = readwstr(f);
    datasets = {};
    fseek(f, dataset_pos, 'bof');
    for i = 1:numdataset
        dataset = readdataset(f);
        fseek(f, dataset.nextdataset, 'bof');
        datasets{end + 1} = dataset;
    end;
    group.datasets = datasets;
return

function dataset = readdataset(f, group)
    datapos = fread(f, 1, 'uint32', 'ieee-be');
    dataset.nextdataset = fread(f, 1, 'uint32', 'ieee-be');
    dataset.name = readwstr(f);

    dataset = readparams(f, dataset);
    dataset.ncols = fread(f, 1, 'uint32', 'ieee-be');

    cols = struct();
    total_typesize = 0;
    for i = 1:dataset.ncols
        cols(i).name = readwstr(f);
        cols(i).value_type = fread(f, 1, 'uint8', 'ieee-be');
        cols(i).type_size = fread(f, 1, 'int32', 'ieee-be');
        total_typesize = total_typesize + cols(i).type_size;
    end;
    dataset.nrows = fread(f, 1, 'uint32', 'ieee-be');
    
    start = ftell(f);

    for i = 1:dataset.ncols
        skip = total_typesize - cols(i).type_size;
        fseek(f, start, 'bof');
        start = start +  cols(i).type_size;
        switch cols(i).value_type
            case 0  %byte
                cols(i).data = fread(f, dataset.nrows, 'int8', skip, 'ieee-be');
            case 1 %ubyte
                cols(i).data = fread(f, dataset.nrows, 'uint8', skip, 'ieee-be');
            case 2 %short
                cols(i).data = fread(f, dataset.nrows, 'int16', skip, 'ieee-be');
            case 3 %ushort
                cols(i).data = fread(f, dataset.nrows, 'uint16', skip, 'ieee-be');
            case 4 %int
                cols(i).data = fread(f, dataset.nrows, 'int32', skip, 'ieee-be');
            case 5 %uint
                cols(i).data = fread(f, dataset.nrows, 'uint32', skip, 'ieee-be');
            case 6 %float
                cols(i).data = fread(f, dataset.nrows, 'float32' ,skip, 'ieee-be');
            case 7 %string
                1/0
            case 8 %wstring
                1/0
        end;
    end;        
    dataset.cols = cols;
return

function s = readstr(f)
    length = fread(f, 1, 'int32', 'ieee-be');
    s = char(fread(f, length, 'char'))';
return

function s = readwstr(f)
    length = fread(f, 1, 'int32', 'ieee-be');
    s = char(fread(f, length, 'uint16', 'ieee-be')');
return

function header = readparams(f, header)
    namevaluenr = fread(f, 1, 'int32', 'ieee-be');
    for i = 1:namevaluenr
        name = strrep(readwstr(f), '-', '_');
        name = strrep(name, '%', '');
            name = strrep(name, 'affymetrix_', '');
        value = uint8(readstr(f));
        type = readwstr(f);
        if(length(type)>100)
            warning('Possible corruption in file, could not find correct type for header parameter');
            header.error_occurred = 1;
            return;
        end;
        switch(type)
            case 'text/x-calvin-integer-8'
                value = swapbytes(typecast(value, 'int8'));
            case 'text/x-calvin-unsigned-integer-8'
                value = swapbytes(typecast(value, 'uint8'));
            case 'text/x-calvin-integer-16'
                value = swapbytes(typecast(value, 'int16'));
            case 'text/x-calvin-unsigned-integer-16'
                value = swapbytes(typecast(value, 'uint16'));
            case 'text/x-calvin-integer-32'
                value = swapbytes(typecast(value, 'int32'));
            case 'text/x-calvin-unsigned-integer-32'
                value = swapbytes(typecast(value, 'uint32'));
            case 'text/x-calvin-float'
                value = swapbytes(typecast(value, 'single'));
            case 'text/plain'
                value = char(swapbytes(typecast(value, 'uint16')));
            case 'text/ascii'
                value = char(value);
        end;
        %header
        %namevaluenr
        %i
        header = setfield(header, name, value);
        %name
        %value
        %type
        %uint16(type)

    end;


function [row,cols]=enterdim(itemnb)
    cols=0;
    rows=0;
    while cols*rows~=itemb
        dimensions=inputdlg({'row nb','col nb'},'enter the chip dimensions',1,{'1164','1164'});
        if ~isempty(dimensions)
            rows=str2num(dimensions{1});
            cols=str2num(dimensions{2});
            if cols*rows~=itemb
                h=warndlg(sprintf('calculated row*cols (%u*%u=%u) ~= %u',cols,cols,cols*cols,itemb));
                waitfor(h)
            end
        else
            h=errodlg('process canceled');
            waitfor(h)
            error('process canceled')
        end
    end
