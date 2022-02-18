function t = rmprefix(s, prefix)
%RMPREFIX removes the prefix given in `prefix` from the string `s` and then return the resulting
% string. It returns `s` if `s` does not end with `prefix`.

t = regexprep(s, ['^', prefix], '');

return
