import numpy as np
from prima import NonlinearConstraint, process_nl_constraints, minimize
import pytest


@pytest.mark.parametrize("lb1", (-np.inf, [-np.inf], np.array([-np.inf])))
@pytest.mark.parametrize("lb2", (-np.inf, [-np.inf]*2, np.array([-np.inf]*2)))
@pytest.mark.parametrize("ub1", (0, [0], np.array([0])))
@pytest.mark.parametrize("ub2", (0, [0]*2, np.array([0]*2)))
def test_multiple_nl_constraints_various_data_types(lb1, ub1, lb2, ub2):
    nlc1 = NonlinearConstraint(lambda x: x, lb=lb1, ub=ub1)
    nlc2 = NonlinearConstraint(lambda x: [x, x], lb=lb2, ub=ub2)
    nlcs = [nlc1, nlc2]
    x0 = 1
    constr_func, nlconstr0 = process_nl_constraints(nlcs, x0)
    assert all(nlconstr0 == [x0, x0, x0])
    assert len(constr_func(0)) == 3
    assert all(constr_func(0) == [0, 0, 0])


@pytest.mark.parametrize("lb", (-np.inf, [-np.inf]*2, np.array([-np.inf]*2)))
@pytest.mark.parametrize("ub", (0, [0]*2, np.array([0]*2)))
def test_single_nl_constraint(lb, ub):
    nlc = NonlinearConstraint(lambda x: [x, x], lb=lb, ub=ub)
    x0 = 2.1
    constr_func, nlconstr0 = process_nl_constraints([nlc], x0)
    assert all(nlconstr0 == [x0, x0])
    assert len(constr_func(0)) == 2
    assert all(constr_func(x0) == [x0, x0])

@pytest.mark.parametrize("lb", (-np.inf, [-np.inf]))
@pytest.mark.parametrize("ub", (0.0, [0.0]))
def test_length_nlc_fun_not_equal_to_length_lb_ub(lb, ub):
    if lb == -np.inf and ub == 0.0:
        return  # No exception here
    nlc = NonlinearConstraint(lambda x: [x, x], lb=lb, ub=ub)
    x0 = 2.1
    with pytest.raises(ValueError) as e_info:
        process_nl_constraints([nlc], x0)
    if lb == -np.inf:
        assert e_info.match("The number of elements in the constraint function's output does not match the number of elements in the upper bound.")
    else:
        assert e_info.match("The number of elements in the constraint function's output does not match the number of elements in the lower bound.")


def test_length_nlc_fun_ne_to_length_ub():
    # This specific test is for the case that both lb and ub are vectors, and lb has the correct
    # length but ub does not.
    nlc = NonlinearConstraint(lambda x: [x, x], lb=[-np.inf, 0], ub=[0])
    x0 = 2.1
    with pytest.raises(ValueError) as e_info:
        process_nl_constraints([nlc], x0)
    assert e_info.match("The number of elements in the constraint function's output does not match the number of elements in the upper bound.")


def test_lb_neg_inf_ub_inf_raises():
    nlc = NonlinearConstraint(lambda x: [x, x], lb=-np.inf, ub=np.inf)
    with pytest.raises(ValueError) as e_info:
        process_nl_constraints([nlc], 0)
    assert e_info.match("A NonlinearConstraint was provided without specifying lower or upper bounds")

def test_lb_neg_inf_ub_vector_w_inf_raises():
    nlc = NonlinearConstraint(lambda x: [x, x], lb=-np.inf, ub=[0, np.inf])
    with pytest.raises(ValueError) as e_info:
        process_nl_constraints([nlc], 0)
    assert e_info.match("A NonlinearConstraint was provided without specifying lower or upper bounds")

def test_lb_vector_with_neg_inf_ub_inf_raises():
    nlc = NonlinearConstraint(lambda x: [x, x], lb=[-np.inf, 0], ub=np.inf)
    with pytest.raises(ValueError) as e_info:
        process_nl_constraints([nlc], 0)
    assert e_info.match("A NonlinearConstraint was provided without specifying lower or upper bounds")

def test_lb_vector_with_neg_inf_ub_vector_w_inf_at_same_index_raises():
    nlc = NonlinearConstraint(lambda x: [x, x], lb=[-np.inf, 0], ub=[np.inf, 1])
    with pytest.raises(ValueError) as e_info:
        process_nl_constraints([nlc], 0)
    assert e_info.match("A NonlinearConstraint was provided without specifying lower or upper bounds")
