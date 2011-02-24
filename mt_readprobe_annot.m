%MT_READPROBE_ANNOT - Read an Affymetrix Probe Annotation file
%
%  PA = MT_READPROBE_ANNOT(FNAME)
%
% INPUT
%   FNAME			File name 
%
% OUTPUT
%		PA				Structure containing probe sequences and position
%
% DESCRIPTION
% Reads in an Affymetrix probe anntoation file. All fields are stored 
% in a self-explanatory structure. 
% Assumes ascii or UTF8 format. In case UTF16 was used, please use the 
% iconv utility to convert the file: 
% iconv -f UTF16 -t UTF8 probetab_file > new_probetab_file
%
% SEE ALSO
% MT_READCDF, MT_READCEL, MT_READGIN, MT_CEL2PROBES

% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function probeannot = mt_readprobe_annot(fname)

   f = fopen(fname,'r');
   if (f <= 0), error ('cannot open file'); end;
   line = fgetl(f);
   line == sprintf('\t')
   nfields = sum(line == sprintf('\t')) + 1;
   frewind(f);

   fdescriptor = '%s%d%d%d%s%s';
   nfields = nfields - 6;
   while(nfields > 0)
      fdescriptor = strcat(fdescriptor,'%s');
      nfields = nfields - 1;
   end;
   res = textscan(f,fdescriptor,'headerLines',1,'Delimiter','\t');
   fclose(f);

   probeannot.probeind = zeros(max(res{2}) + 1,max(res{3}) + 1);
   probeannot.name = res{1};
   probeannot.pos = [res{2} + 1 res{3} + 1];
   probeannot.position = res{4};
   probeannot.sequence = res{5};

   fprintf(1,'Mapping probesets\n');
   pos = probeannot.pos;
   for i = 1:length(res{1})
      probeannot.probeind(pos(i,1),pos(i,2)) = i;
   end;

