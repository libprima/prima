function is_int = isint(x)
is_int = (abs(x-round(x)) < eps);
