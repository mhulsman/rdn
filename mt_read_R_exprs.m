%MT_READ_R_EXPRS - Read summarized data from R data file

%  E = MT_WRITE_R_EXPRS(probes,filename)
%
% INPUT
%   PROBES			Probe structure into which to load summarized data
%                   Overwrites overall_factors and array_factors,
%                   matches array names in R file to array names in probes structure,
%                   and probeset names to probeset names in probes structure.
%   FILENAME        file to read from
%
% SEE ALSO
% MT_WRITE_R_EXPRS
%
% (c) Marc Hulsman, 2010
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands
function probes = mt_read_R_exprs(probes,fname)

%read array names
f = fopen(fname,'r');
frewind(f);
tmp = fgetl(f);
%process array names
start = regexp(tmp,'([\w\.]+)');
if(length(strfind(lower(tmp),'.cel')) > 0)
   stop = [start(2:end) - 6 length(tmp) - 4];
else
   stop = [start(2:end) - 2 length(tmp)];
end;
array_names = {};
for i = 1:length(start)
   array_names{i} = tmp(start(i):stop(i));
end;

%read probe names
tmp = textscan(f,'%s%*[^\n]');
probe_names = tmp{1};
fclose(f);

%read data
values = dlmread(fname,'\t',1,1);

%matches probes names
probes_idx = [];
for i = 1:length(probe_names);
   probes_idx(find(strcmp(probe_names{i},probes.name))) = i;
end;

%matches array names
acounter = 1;
for i = 1:length(array_names)
   t = find(strcmp(array_names{i},probes.array_names));
   if(length(t) > 0)
      array_idx(t) = acounter;
      acounter = acounter + 1;
   end;
end;

%determine data
sv = values(probes_idx,array_idx)';
%perform log2 scaling if necessary
if(max(max(sv)) > 50)
   sv = log2(sv);
end;

%set variables
probes.overall_factors = mean(sv);
probes.array_factors = sv - repmat(probes.overall_factors,length(probes.array_names),1);
