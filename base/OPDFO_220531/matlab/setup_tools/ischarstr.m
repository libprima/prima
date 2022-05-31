function iscs = ischarstr(x)
%ISCHARSTR checks whether an input is a `char` or `string`

iscs = (isa(x, 'char') || isa(x, 'string'));
