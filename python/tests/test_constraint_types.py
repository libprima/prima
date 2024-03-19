from prima import minimize, LinearConstraint as LC, NonlinearConstraint as NLC
from scipy.optimize import LinearConstraint as ScipyLC, NonlinearConstraint as ScipyNLC
import numpy as np
from objective import fun
import pytest

def test_prima_lc_is_scipy_lc():
    assert LC is ScipyLC


def test_linear_constraint_object():
    lc = LC(np.array([[1, 1],[1, -1]]), lb=[0, 0], ub=[8, 2])
    x0 = [0, 0]
    res = minimize(fun, x0, constraints=lc)
    assert np.isclose(res.x[0], 4.5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 3.5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 0.5, atol=1e-6, rtol=1e-6)


def test_linear_constraint_dict():
    lc = {'A': np.array([[1, 1],[1, -1]]), 'lb':[0, 0], 'ub':[8, 2]}
    x0 = [0, 0]
    res = minimize(fun, x0, constraints=lc)
    assert np.isclose(res.x[0], 4.5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 3.5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 0.5, atol=1e-6, rtol=1e-6)
    

def test_prima_nlc_is_scipy_nlc():
    assert NLC is ScipyNLC


def test_nonlinear_constraint_object():
    nlc = NLC(lambda x: np.array([x[0]**2, x[1]**2]), lb=[25]*2, ub=[100]*2)
    x0 = [0, 0]
    res = minimize(fun, x0, constraints=nlc)
    assert np.isclose(res.x[0], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 1, atol=1e-6, rtol=1e-6)


def test_nonlinear_constraint_dict():
    nlc = {'fun': lambda x: np.array([x[0]**2, x[1]**2]), 'lb':[25]*2, 'ub':[100]*2}
    x0 = [0, 0]
    res = minimize(fun, x0, constraints=nlc)
    assert np.isclose(res.x[0], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 1, atol=1e-6, rtol=1e-6)


def test_unsupported_type_raises_exception():
    # By not providing ub this becomes an unsupported type
    lc = {'A': np.array([1, 1]), 'lb': 0}
    x0 = [0, 0]
    
    with pytest.raises(ValueError) as e_info:
        minimize(fun, x0, constraints=lc)
    assert e_info.match('Constraint type not recognized')
