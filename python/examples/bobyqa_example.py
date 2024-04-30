from prima import minimize
from objective import fun

x0 = [0.0] * 2

# Unbounded the min is (5, 4), so we put a bound on the first
# element to make sure we're exercising the bounding logic.
bounds = [(None, 2.5), (-6, 6)]  # x1 from -inf to 2.5, x2 from -6 to 6
res = minimize(fun, x0, bounds=bounds)
print(res)
assert abs(res.x[0] - 2.5) < 2e-1 and abs(res.x[1] - 4.0) < 2e-1

# Test not all bounds being provided
bounds = [(None, 2.5)]
res = minimize(fun, x0, method='BOBYQA', bounds=bounds)
print(res)
assert abs(res.x[0] - 2.5) < 2e-1 and abs(res.x[1] - 4.0) < 2e-1

# Test bounding only the second variable
bounds = [None, (None, 6)]
res = minimize(fun, x0, method='BOBYQA', bounds=bounds)
print(res)
assert abs(res.x[0] - 5) < 2e-1 and abs(res.x[1] - 4.0) < 2e-1
