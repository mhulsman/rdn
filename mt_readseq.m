%READSEQ - Read an Affymetrix probeset sequence file
%
%  SEQ = READSEQ(FNAME,SEQ)
%
% INPUT
%   FNAME			File name 
%   SEQ             [Optional] Probeset seq structure containing already
%                              sequences from other probeset
%                              sequence file. 
%
% OUTPUT
%		SEQ				Structure containing probeset sequences from file and 
%                       the optional SEQ parameter
%
% DESCRIPTION
% Reads in an Affymetrix probeset sequence file. All fields are stored 
% in a self-explanatory structure. 
%
% SEE ALSO
% READCDF, READCEL, READGIN, READPROBE_ANNOT, CEL2PROBES

% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function seq = mt_readseq(fname,seq,progress)

	if (nargin < 3), progress = 1; end;

	f = fopen(fname,'r');
	if (f <= 0), error ('cannot open file'); end;
	frewind(f);


   if(nargin < 2 | ~isfield(seq,'name'))
      seq.name = {};
      seq.sequence = {};
   end;

   done = 0; 
   k = length(seq.name) + 1;
   line = fgets(f);

   fprintf(1,'Reading sequence file: ');
	while (~feof(f))
      if(mod(k,1000) == 0) fprintf(1,'.'); end;
      n = regexp(line,'^>\w+\:[\w-_]+\:(?<probeset>[\w-_/\.]+)\;','names');
      sequence = '';
      line = 'x';
      while((~feof(f)) & (line(1) ~= '>'))
         line = fgetl(f);
         while((~feof(f)) & (length(line) == 0))
            line = fgetl(f);
         end;
         if((length(line) > 0) & (line(1) ~= '>'))
            sequence = strcat(sequence,line);
         end;
      end;
      if(length(find(strcmp(n.probeset, seq.name))) == 0) 
          seq.name{k} = n.probeset;
          seq.sequence{k} = sequence;
          k = k + 1;
      end;
	end;
   fprintf(1,'\n');

   fclose(f);

return

