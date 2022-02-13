function iss = is_str(x)
%IS_STR checks whether an input is a `char` or `string`

iss = (isa(x, 'char') || isa(x, 'string'));

return
