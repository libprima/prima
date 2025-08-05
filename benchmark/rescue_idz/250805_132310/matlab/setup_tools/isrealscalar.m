function isrs = isrealscalar(x)
%ISREALSCALAR checks whether x is a real scalar.
% N.B.: isrealscalar([]) = FALSE, isrealscalar(NaN) = TRUE, isrealscalar(inf) = TRUE!!!

isrs = isnumeric(x) && isreal(x) && isscalar(x);

return
