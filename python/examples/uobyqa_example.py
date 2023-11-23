from prima import minimize
from objective import fun


x0 = [0.0] * 2


res = minimize(fun, x0, method='UOBYQA')
print(res)
assert fun.result_point_and_value_are_optimal(res)
