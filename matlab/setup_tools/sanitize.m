function sanitize(directory)
%SANITIZE cleans up `directory` so that it is proper for the compilation. Without doing this, files
% may be linked mistakenly, leading to runtime errors such as SEGFAULT.

cellfun(@(filename) delete(filename), list_modo_files(directory));

return
