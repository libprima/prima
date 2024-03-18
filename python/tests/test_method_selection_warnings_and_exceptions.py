from prima import minimize, LinearConstraint, NonlinearConstraint, Bounds
import pytest


def test_method_not_recognized():
    with pytest.raises(ValueError) as e_info:
        minimize(lambda x: x, [0.0], method="not a method")
    assert e_info.match("Method must be one of NEWUOA, UOBYQA, BOBYQA, COBYLA, or LINCOA, not 'not a method'")
    

@pytest.mark.parametrize("method", ["newuoa", "uobyqa", "bobyqa", "lincoa"])
def test_providing_nonlinear_constraints_to_non_cobyla_method(method):
    with pytest.raises(ValueError) as e_info:
        minimize(lambda x: x, [0.0], method=method, constraints=NonlinearConstraint(lambda x: x, 0, 1))
    assert e_info.match("Nonlinear constraints were provided for an algorithm that cannot handle them")


@ pytest.mark.parametrize("method", ["newuoa", "uobyqa", "bobyqa"])
def test_providing_linear_constraints_to_non_cobyla_non_lincoa_method(method):
    with pytest.raises(ValueError) as e_info:
        minimize(lambda x: x, [0.0], method=method, constraints=LinearConstraint([1], 0, 1))
    assert e_info.match("Linear constraints were provided for an algorithm that cannot handle them")


@pytest.mark.parametrize("method", ["newuoa", "uobyqa"])
def test_providing_bounds_to_non_cobyla_non_lincoa_non_bobyqa_method(method):
    with pytest.raises(ValueError) as e_info:
        minimize(lambda x: x, [0.0], method=method, bounds=Bounds([0], [1]))
    assert e_info.match("Bounds were provided for an algorithm that cannot handle them")
