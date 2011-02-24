%MT_READGENE_ANNOT - Read an Affymetrix probeset annotation file
%
%  PROBES = MT_READCDF(PROBES, filename)
%
% INPUT
%   PROBES          Probes structure
%   FNAME			File name 
%
% OUTPUT
%   PROBES			Probes structure with probeset annotations
%
% DESCRIPTION
% Reads in an Affymetrix probeset annotation file. All fields are 
% added to the probes structure.
%
%
% SEE ALSO
% MT_READCEL,MT_READPROBE_ANNOT,MT_READSEQ,MT_READCDF,MT_CEL2PROBES

% (c) Marc Hulsman, 2008
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands
function p = mt_readgene_annot(p,filename);

f = fopen(filename, 'r');
k = textscan(f,'%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%*[^\n]','bufsize',1280000,'delimiter',',','headerlines',1);


names = p.name;
anames = k{1};

s_names = [];
as_names = [];
fprintf(1,'Mapping names');
for i=1:length(names)
    if(mod(i,1000) == 0)
        fprintf(1,'.');
    end;
    res = find(strcmp(anames,names{i}));
    if(length(res) == 0)
        warning('Could not find annotation for probeset %s',names{i});
    else
        as_names(end + 1) = res(1);
        s_names(end + 1) = i;
    end;
end;
fprintf(1,'\n');

p.gene_names = {};
p.gene_description = {};

t = k{15};
p.gene_names(s_names) = t(as_names); 
t = k{14};
p.gene_description(s_names) = t(as_names); 
t = k{13};
p.chrom_loc_det(s_names) = t(as_names);
t = k{16};
p.chrom_loc(s_names) = t(as_names);
fclose(f);
