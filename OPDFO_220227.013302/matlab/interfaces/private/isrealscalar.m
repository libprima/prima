function isrs = isrealscalar(x)
% isrealscalar([]) = FALSE, isrealscalar(NaN) = TRUE, isrealscalar(inf) = TRUE!!!
isrs = isnumeric(x) && isreal(x) && isscalar(x);
return
