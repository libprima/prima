import numpy as np

def fun(x):
    x1, x2 = x
    f = (x1 - 5) ** 2 + (x2 - 4.0) ** 2
    return f

fun.optimal_x = np.array((5, 4))
fun.result_point_is_optimal = lambda result: np.isclose(0, np.linalg.norm(result.x - fun.optimal_x), atol=1e-6, rtol=1e-6)
fun.optimal_f = fun(fun.optimal_x)
fun.result_point_and_value_are_optimal = lambda result: fun.result_point_is_optimal(result) and np.isclose(result.fun, fun.optimal_f, atol=1e-6, rtol=1e-6)
