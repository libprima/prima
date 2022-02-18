function t = rmsuffix(s, suffix)
%RMSUFFIX removes the suffix given in `suffix` from the string `s` and then return the resulting
% string. It returns `s` if `s` does not end with `suffix`.

t = regexprep(s, [suffix, '$'], '');

return
