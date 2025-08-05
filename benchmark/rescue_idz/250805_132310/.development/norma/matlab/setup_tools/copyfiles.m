function copyfiles(file_list, directory)
%COPYFILES copies the files given in `file_list` to `directory`.

cellfun(@(filename) copyfile(filename, directory), file_list);
return
