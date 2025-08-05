function [isrm, m, n] = isrealmatrix(x)
%ISREALMATRIX checks whether x is a real matrix.
% N.B.: isrealmatrix([]) = true

if isempty(x)
    isrm = true;
    m = 0;
    n = 0;
elseif isnumeric(x) && isreal(x) && ismatrix(x)
    isrm = true;
    [m, n] = size(x);
else
    isrm = false;
    m = NaN;
    n = NaN;
end
return
