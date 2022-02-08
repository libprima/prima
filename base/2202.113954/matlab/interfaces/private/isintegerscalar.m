function isis = isintegerscalar(x)
% isintegerscalar([]) = FALSE, isintegerscalar(NaN) = FALSE, isintegerscalar(inf) = FALSE !!!
isis = isrealscalar(x) && (rem(x,1) == 0);
return
