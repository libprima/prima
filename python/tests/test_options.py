from prima import minimize, NonlinearConstraint as NLC, PRIMAMessage
from objective import fun
import numpy as np
import pytest
from sys import platform


def fun_with_star_args(x, *args):
    return fun(x) + args[0][0]


def fun_with_regular_args(x, args):
    return fun(x) + args[0]


def test_normal_function():
    x0 = [0.0] * 2
    res = minimize(fun, x0)
    assert fun.result_point_and_value_are_optimal(res)


def test_function_with_star_args():
    x0 = np.array([0.0] * 2)
    myargs = (5,)
    res = minimize(fun_with_star_args, x0, args=myargs)
    assert fun.result_point_is_optimal(res)
    assert np.isclose(res.fun, fun.optimal_f + myargs[0], rtol=1e-6)


def test_function_with_regular_args():
    x0 = np.array([0.0] * 2)
    myargs = (6,)
    res = minimize(fun_with_regular_args, x0, args=myargs)
    assert fun.result_point_is_optimal(res)
    assert np.isclose(res.fun, fun.optimal_f + myargs[0], rtol=1e-6)


def test_callback():
    x0 = [0.0] * 2
    nlc = NLC(lambda x: np.array([x[0], x[1]]), lb=[-np.inf]*2, ub=[10]*2)
    callback_called = False
    def callback(x, f, nf, tr, cstrv, nlconstr):
        nonlocal callback_called  # this makes sure we reference the variable in the parent scope
        callback_called = True
    res = minimize(fun, x0, method="COBYLA", constraints=nlc, callback=callback)
    assert callback_called
    assert fun.result_point_and_value_are_optimal(res)


def test_callback_early_termination():
    x0 = [0.0] * 2
    def callback(x, *args):
        if x[0] > 1:
            return True
    res = minimize(fun, x0, callback=callback)
    # Assert that the results are not optimal since we terminated early
    assert not fun.result_point_and_value_are_optimal(res)


def test_ftarget():
    x0 = [0.0] * 2
    options = {'ftarget': 20}
    res = minimize(fun, x0, options=options)
    # Assert that the results are not optimal since we're stopping
    # at ftarget=20 instead of the natural minimum of 0.
    assert not fun.result_point_and_value_are_optimal(res)


@pytest.mark.skipif(platform == "win32", reason="Windows outputs some strange characters, probably \\r\\n")
def test_iprint(capfd):
    x0 = [0.0] * 2
    options = {'iprint': PRIMAMessage.EXIT}
    res = minimize(fun, x0, options=options)
    assert fun.result_point_and_value_are_optimal(res)
    outerr = capfd.readouterr()
    assert outerr.out == '''No bounds or constraints detected, applying NEWUOA

Return from NEWUOA because the trust region radius reaches its lower bound.
Number of function values = 23   Least value of F =  0.000000000000000E+000
The corresponding X is:  5.000000000000000E+000   4.000000000000000E+000
'''
    assert outerr.err == ''


def test_maxfun():
    x0 = [0.0] * 2
    options = {'maxfun': 10}
    res = minimize(fun, x0, options=options)
    assert res.nfev == 10


def test_maxfev():
    x0 = [0.0] * 2
    options = {'maxfev': 10}
    res = minimize(fun, x0, options=options)
    assert res.nfev == 10


def test_npt():
    x0 = [0.0] * 2
    options = {'npt': 4}
    res_with_npt = minimize(fun, x0, options=options)
    res_without_npt = minimize(fun, x0)
    assert res_with_npt.nfev != res_without_npt.nfev


def test_rhobeg():
    x0 = [0.0] * 2
    options = {'rhobeg': 1e-3}
    res_with_rhobeg = minimize(fun, x0, options=options)
    res_without_rhobeg = minimize(fun, x0)
    # With a smaller rhobeg, we're taking smaller steps at first
    # and so we should converge more slowly
    assert res_with_rhobeg.nfev > res_without_rhobeg.nfev


def test_rhoend():
    x0 = [0.0] * 2
    options = {'rhoend': 1e-2}
    res_with_rhoend = minimize(fun, x0, options=options)
    res_without_rhoend = minimize(fun, x0)
    # With a smaller rhoend we have a larger tolerance for the final
    # trust region radius, and so we should converge more quickly
    assert res_with_rhoend.nfev < res_without_rhoend.nfev
