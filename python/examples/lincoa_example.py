from prima import minimize, LinearConstraint as LC
from objective import fun
import numpy as np


x0 = [0.0] * 2

# x1 from -inf to 2.5, x2 from -6 to 6
bounds = [(None, 2.5), (-6, 6)]

Aineq = np.array([[1.0, 0.0],
                  [0.0, 1.0],
                  [1.0, 1.0]])
bineq = np.array([2.7, # The bounds is below this, so it should take precedence
                  1.5, # This is below the bounds, so it should take precedence
                  10.0])

constraint = LC(Aineq, -np.inf, bineq)

res = minimize(fun, x0, method='LINCOA', bounds=bounds, constraints=constraint)
print(res)
assert abs(res.x[0] - 2.5) < 2e-1 and abs(res.x[1] - 1.5) < 2e-1
