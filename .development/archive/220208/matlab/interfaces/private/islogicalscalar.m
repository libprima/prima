function isls = islogicalscalar(x) % islogicalscalar([]) = FALSE !!!
if isa(x, 'logical') && isscalar(x)
    isls = true;
elseif isrealscalar(x) && (x==1 || x==0) % !!!!!!
    isls = true;
else
    isls = false;
end
return
