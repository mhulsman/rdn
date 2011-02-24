% get all filenames in 'directory', adhering to regular expression expr_lower  
% (applied to lowercase version of the filename)
% use filenames = dir_filter(dirpath, '.cel') to get all cell files in a directory.

function file_names = dir_filter(directory,expr_lower)

x = dir(directory);

file_names = {};
k = 0;

for i = 1:length(x)
   h = regexp(lower(x(i).name),expr_lower);
   if(~isempty(h))
      k = k + 1;
      file_names{k} = strcat(directory,x(i).name);
   end;
end;

