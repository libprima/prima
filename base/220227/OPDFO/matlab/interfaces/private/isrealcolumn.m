function [isrc, len] = isrealcolumn(x) % isrealcolumn([]) = true
if isempty(x)
    isrc = true;
    len = 0;
elseif isnumeric(x) && isreal(x) && isvector(x) && (size(x, 2) == 1)
    isrc = true;
    len = length(x);
else
    isrc = false;
    len = NaN;
end
return
