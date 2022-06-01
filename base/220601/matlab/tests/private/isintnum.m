function is_intnum = isintnum(x)
is_intnum = isnumeric(x) && numel(x) == 1 && abs(x-round(x)) < eps;
