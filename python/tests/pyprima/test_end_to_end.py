from prima import minimize, Bounds, LinearConstraint
import numpy as np

def obj(x):
    return (x[0] - 1)**2 + (x[1] - 2.5)**2
obj.x0 = np.array([5, 5])
obj.optimal = np.array([1, 2.5])


def test_end_to_end_bounds(backend_fixture):
    bounds = Bounds([-5, 10], [5, 10])
    result = minimize(obj, obj.x0, method='cobyla', bounds=bounds, options={'backend': backend_fixture})
    assert np.allclose(result.x, np.array([1, 10]), atol=1e-3)


def test_end_to_end_linear_constraints(minimize_with_debugging, backend_fixture):
    # x1 + x2 = 5
    A = np.array([[1, 1]])
    b = np.array([5])
    result = minimize_with_debugging(obj, obj.x0, method='cobyla', constraints=[LinearConstraint(A, b, b)], options={'backend': backend_fixture})
    assert np.allclose(result.x, np.array([1.75, 3.25]), atol=1e-3)
    assert np.allclose(A @ result.x, b, atol=1e-3)
