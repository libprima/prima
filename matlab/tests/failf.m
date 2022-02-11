function f = failf(x)

f = sin(x);

if (abs(x) > 1)
    error('Fail at %f\n', x);
end
