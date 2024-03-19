from prima import minimize, NonlinearConstraint as NLC
from objective import fun
import numpy as np  # To set lb to -inf


def f_con(x):
    con = x[0]**2 - 9
    return con


x0 = [0.0] * 2
myNLC = NLC(f_con, [-np.inf], [0])

res = minimize(fun, x0, method='COBYLA', constraints=myNLC)
print(res)
assert abs(res.x[0] - 3) < 1e-2 and abs(res.x[1] - 4) < 1e-2 and abs(res.fun - 4) < 1e-2
