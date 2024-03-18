from prima import minimize, NonlinearConstraint as NLC, LinearConstraint as LC
from objective import fun
import numpy as np
import pytest

def test_x0_as_list():
    x0 = [0.0] * 2
    res = minimize(fun, x0)
    assert fun.result_point_and_value_are_optimal(res)


def test_x0_as_array():
    x0 = np.array([0.0] * 2)
    res = minimize(fun, x0)
    assert fun.result_point_and_value_are_optimal(res)


def test_x0_as_scalar():
    x0 = 0.0
    # We need a custom function since the default objective function we're
    # using for tests expects something that can be unpacked into two variables.
    res = minimize(lambda x: (x-5)**2, x0)
    assert np.isclose(res.x, 5.0, rtol=1e-6)
    assert np.isclose(res.fun, 0.0, rtol=1e-6)


def test_constraint_function_returns_numpy_array():
    nlc = NLC(lambda x: np.array([x[0], x[1]]), lb=[-np.inf]*2, ub=[10]*2)
    x0 = [0, 0]
    res = minimize(fun, x0, method='COBYLA', constraints=nlc)
    assert fun.result_point_and_value_are_optimal(res)


def test_constraint_function_returns_list():
    nlc = NLC(lambda x: [x[0], x[1]], lb=[-np.inf]*2, ub=[10]*2)
    x0 = [0, 0]
    res = minimize(fun, x0, method='COBYLA', constraints=nlc)
    assert fun.result_point_and_value_are_optimal(res)


def test_constraint_function_returns_scalar():
    nlc = NLC(lambda x: float(np.linalg.norm(x) - np.linalg.norm(fun.optimal_x)), lb=[-np.inf], ub=[0])
    x0 = [0, 0]
    res = minimize(fun, x0, method='COBYLA', constraints=nlc)
    assert fun.result_point_and_value_are_optimal(res)



@pytest.mark.parametrize("A", (1, [1], np.array([1])))
@pytest.mark.parametrize("lb", (0, [0], np.array([0])))
@pytest.mark.parametrize("ub", (4, [4], np.array([4])))
def test_linear_constraint_data_types(A, lb, ub):
    myLC = LC(A=A, lb=lb, ub=ub)
    # If A is scalar, x must have dimension 1, so we need a univariate function for that
    scalar_fun = lambda x: (x - 5)**2
    x0 = 0
    res = minimize(scalar_fun, x0, method='LINCOA', constraints=myLC)
    assert np.isclose(res.x, 4, rtol=1e-2)


def test_nonlinear_constraint_lb_scalar_ub_scalar():
    nlc = NLC(lambda x: [x[0]**2, x[1]**2], lb=25, ub=100)
    x0 = [0, 0]
    res = minimize(fun, x0, constraints=nlc)
    assert np.isclose(res.x[0], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 1, atol=1e-6, rtol=1e-6)


def test_nonlinear_constraint_lb_scalar_ub_not_scalar():
    nlc = NLC(lambda x: [x[0]**2, x[1]**2], lb=9, ub=[9, 16])
    x0 = [0, 0]
    res = minimize(fun, x0, constraints=nlc)
    assert np.isclose(res.x[0], 3, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 4, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 4, atol=1e-6, rtol=1e-6)


def test_nonlinear_constraint_lb_not_scalar_ub_scalar():
    nlc = NLC(lambda x: [x[0]**2, x[1]**2], lb=[25, 36], ub=100)
    x0 = [0, 0]
    res = minimize(fun, x0, constraints=nlc)
    assert np.isclose(res.x[0], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 6, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 4, atol=1e-6, rtol=1e-6)
    

def test_nonlinear_constraint_lb_not_scalar_ub_not_scalar():
    nlc = NLC(lambda x: [x[0]**2, x[1]**2], lb=[4, 4], ub=[4, 9])
    x0 = [0, 0]
    res = minimize(fun, x0, constraints=nlc)
    assert np.isclose(res.x[0], 2, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 3, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 10, atol=1e-6, rtol=1e-6)
