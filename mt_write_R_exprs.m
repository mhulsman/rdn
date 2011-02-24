%MT_WRITE_R_EXPRS - Write summarized expression to format suitable for R
%
%  E = MT_WRITE_R_EXPRS(probes,filename)
%
% INPUT
%   PROBES			Probe structure 
%   FILENAME        file to write to
%
% SEE ALSO
% MT_READ_R_EXPRS
%
% (c) Marc Hulsman, 2010
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands
function mt_write_R_exprs(probes,filename)

%write first row with array names
f = fopen(filename,'w');
arraynames = probes.array_names;
row1 = sprintf(repmat('%s,',1,size(arraynames,1)),arraynames{:});
fwrite(f,sprintf(','));
fwrite(f,row1(1:(end-1)));
fwrite(f,sprintf('\n'));

%calculate signal
signal = repmat(probes.overall_factors,size(probes.array_factors,1),1) + probes.array_factors;

%write data for each probe
for i = 1:length(probes.name)
    fwrite(f,[char(probes.name{i}),sprintf(',')]);
    col = sprintf('%f,',signal(:,i));
    fwrite(f,col(1:(end-1)));
    fwrite(f,sprintf('\n'));
end;


fclose(f)

