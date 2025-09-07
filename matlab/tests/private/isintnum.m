function is_intnum = isintnum(x)
is_intnum = isnumeric(x) && isscalar(x) && abs(x-round(x)) <= eps;
