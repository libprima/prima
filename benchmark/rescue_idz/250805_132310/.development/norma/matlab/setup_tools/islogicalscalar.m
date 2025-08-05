function isls = islogicalscalar(x)
%ISLOGICALSCALAR checks whether x is a logical scalar, including 0 and 1.
% N.B.: islogicalscalar([]) = FALSE !!!

if isa(x, 'logical') && isscalar(x)
    isls = true;
elseif isrealscalar(x) && (x == 1 || x == 0) % !!!!!!
    isls = true;
else
    isls = false;
end

return
