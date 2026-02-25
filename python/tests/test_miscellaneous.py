from prima import minimize, NonlinearConstraint
from prima.backends.pyprima.common.infos import CALLBACK_TERMINATE, SMALL_TR_RADIUS
from prima.backends.pybindings import __cmake_build_type__
import numpy as np
import pytest
import platform

def obj(x):
    return (x[0] - 1)**2 + (x[1] - 2.5)**2
obj.x0 = np.array([5, 5])
obj.optimal = np.array([1, 2.5])


def test_callback_terminate(pyprima_turn_on_debugging, backend_fixture):
    def callback(x, *args):
        return True
    result = minimize(obj, obj.x0, method='cobyla', callback=callback, options={'backend': backend_fixture})
    assert result.nfev == 3
    assert result.status == CALLBACK_TERMINATE


def test_callback_no_terminate(backend_fixture):
    def callback(x, *args):
        pass
    result = minimize(obj, obj.x0, method='cobyla', callback=callback, options={'backend': backend_fixture})
    # Different platforms finish with different nfev due to different optimizations
    assert result.nfev == (56 if (
        (platform.machine().lower() in ["x86_64", "amd64", "i386"]) or
        (__cmake_build_type__ == "Debug") or (backend_fixture == 'Python')
        ) else 54)
    assert np.allclose(result.x, obj.optimal, atol=1e-3)
    assert result.status == SMALL_TR_RADIUS


def test_rhoend_without_rhobeg(backend_fixture):
    result = minimize(obj, obj.x0, method='cobyla', options={'rhoend': 4e-4, 'backend': backend_fixture})
    assert np.allclose(result.x, obj.optimal, atol=1e-3)


def test_rhobeg_without_rhoend(backend_fixture):
    # Needs to be negative to trigger the right section of code.
    if backend_fixture == "Python":
        with pytest.warns(UserWarning, match="COBYLA: Invalid RHOBEG; it should be a positive number; it is set to 1"):
            result = minimize(obj, obj.x0, method='cobyla', options={'rhobeg': -1, 'backend': backend_fixture})
    else:
        # The Fortran emits warnings to stdout as opposed to a Python warning
        result = minimize(obj, obj.x0, method='cobyla', options={'rhobeg': -1, 'backend': backend_fixture})
    assert np.allclose(result.x, obj.optimal, atol=1e-3)


def test_eta2_without_eta1(backend_fixture):
    result = minimize(obj, obj.x0, method='cobyla', options={'eta2': 0.7, 'backemd': backend_fixture})
    assert np.allclose(result.x, obj.optimal, atol=1e-3)


def test_eta2_without_eta1_and_eta2_out_of_range(backend_fixture):
    if backend_fixture == 'Python':
        with pytest.warns(UserWarning, match=r"COBYLA: Invalid ETA2; it should be in the interval \[0, 1\) and not less than ETA1; it is set to 0.7000000000000001"):
            result = minimize(obj, obj.x0, method='cobyla', options={'eta2': 1.7, 'backend': backend_fixture})
    else:
        result = minimize(obj, obj.x0, method='cobyla', options={'eta2': 1.7, 'backend': backend_fixture})
    assert np.allclose(result.x, obj.optimal, atol=1e-3)


def test_eta1_without_eta2_and_eta1_out_of_range(backend_fixture):
    if backend_fixture == 'Python':
        with pytest.warns(UserWarning, match=r"COBYLA: Invalid ETA1; it should be in the interval \[0, 1\) and not more than ETA2; it is set to 0.09999999999999999"):
            with pytest.warns(UserWarning, match=r"COBYLA: Invalid ETA2; it should be in the interval \[0, 1\) and not less than ETA1; it is set to 0.7000000000000001"):
                result = minimize(obj, obj.x0, method='cobyla', options={'eta1': 1.7, 'backend': backend_fixture})
    else:
        result = minimize(obj, obj.x0, method='cobyla', options={'eta1': 1.7, 'backend': backend_fixture})
    assert np.allclose(result.x, obj.optimal, atol=1e-3)

@pytest.mark.parametrize('iprint', [-1, 0, 1, 2, 3, 4])
def test_iprint(iprint, backend_fixture):
    if backend_fixture == 'Python' and iprint == 4:
        with pytest.warns(UserWarning, match=r"COBYLA: Invalid IPRINT; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to 0"):
            result = minimize(obj, obj.x0, method='cobyla', options={'iprint': iprint, 'backend': backend_fixture})
    else:
        result = minimize(obj, obj.x0, method='cobyla', options={'iprint': iprint, 'backend': backend_fixture})
    assert np.allclose(result.x, obj.optimal, atol=1e-3)


def test_minimize_constraint_violation(backend_fixture):
    # We set up conflicting constraints so that the algorithm will be
    # guaranteed to end up with maxcv > 0.
    cons = [NonlinearConstraint(lambda x: x - 4, -np.inf, 0),
            NonlinearConstraint(lambda x: 5 - x, -np.inf, 0)]
    result = minimize(lambda x: x[0], np.array([0]), method='cobyla', constraints=cons,
                    options={'backend': backend_fixture})
    assert result.maxcv > 0.1
    assert result.status == SMALL_TR_RADIUS


def test_scalar():
    result = minimize(lambda x: x**2, 5, method='cobyla')
    assert np.allclose(result.x, 0, atol=1e-3)