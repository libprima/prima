from prima import minimize
from objective import fun

def callback(x, f, nf, tr, *args):
    print(x, f, nf, tr)
    return bool(x[0] > 3)  # Testing terminate functionality

x0 = [0.0] * 2

# Test default arguments
res = minimize(fun, x0, method='NEWUOA')
print(res)  # test repr
assert fun.result_point_and_value_are_optimal(res)

# Test callback and options
options = {'rhobeg': 0.1}
res = minimize(fun, x0, method='NEWUOA', callback=callback, options=options)
print(res.message)
assert not fun.result_point_and_value_are_optimal(res)
