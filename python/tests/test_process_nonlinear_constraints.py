import numpy as np
from prima import NonlinearConstraint, process_nl_constraints
import pytest


@pytest.mark.parametrize("lb1", (-np.inf, [-np.inf], np.array([-np.inf])))
@pytest.mark.parametrize("lb2", (-np.inf, [-np.inf]*2, np.array([-np.inf]*2)))
@pytest.mark.parametrize("ub1", (0, [0], np.array([0])))
@pytest.mark.parametrize("ub2", (0, [0]*2, np.array([0]*2)))
def test_multiple_nl_constraints_various_data_types(lb1, ub1, lb2, ub2):
    nlc1 = NonlinearConstraint(lambda x: x, lb=lb1, ub=ub1)
    nlc2 = NonlinearConstraint(lambda x: [x, x], lb=lb2, ub=ub2)
    nlcs = [nlc1, nlc2]
    constr_func = process_nl_constraints(nlcs)
    assert all(constr_func(0) == [0, 0, 0])


@pytest.mark.parametrize("lb", (-np.inf, [-np.inf], np.array([-np.inf])))
@pytest.mark.parametrize("ub", (0, [0], np.array([0])))
def test_single_nl_constraint(lb, ub):
    nlc = NonlinearConstraint(lambda x: [x, x], lb=lb, ub=ub)
    x0 = 2.1
    constr_func = process_nl_constraints([nlc])
    assert all(constr_func(x0) == [x0, x0])
