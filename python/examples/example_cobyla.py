import prima
from scipy.optimize import Bounds

def fun(x):
    x1, x2 = x
    f = 5.0 * (x1 - 3.0) ** 2 + 7.0 * (x2 - 2.0) ** 2 + 0.1 * (x1 + x2) - 10.0
    return f


def f_con(x):
    x1, x2 = x
    # ||x||^2<=13
    con = [x1**2 + x2**2 - 13]
    return con


x0 = [0.0] * 2
xl = [-6.0] * 2
xu = [6.0] * 2
bounds = Bounds(xl, xu)
res = prima.cobyla(fun, x0, f_con, 1, bounds=bounds)
print(res)
assert abs(res.x[0] - 3.0) < 2e-1 and abs(res.x[1] - 2.0) < 2e-1