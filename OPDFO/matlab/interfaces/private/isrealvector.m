function [isrv, len] = isrealvector(x)  % isrealvector([]) = true
if isrealrow(x) || isrealcolumn(x)
    isrv = true;
    len = length(x);
else
    isrv = false;
    len = NaN;
end
return
