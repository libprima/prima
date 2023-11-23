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
    x0 = 0
    options = {}
    nlc = process_nl_constraints(x0, nlcs, options)
    assert all(nlc.lb == [-np.inf, -np.inf, -np.inf])
    assert all(nlc.ub == [0, 0, 0])
    assert all(nlc.fun(0) == [0, 0, 0])
    assert options['nlconstr0'] == [0, 0, 0]


@pytest.mark.parametrize("lb", (-np.inf, [-np.inf], np.array([-np.inf])))
@pytest.mark.parametrize("ub", (0, [0], np.array([0])))
def test_single_nl_constraint(lb, ub):
    nlc = NonlinearConstraint(lambda x: [x, x], lb=lb, ub=ub)
    options = {}
    x0 = 2.1
    processed_nlc = process_nl_constraints(x0, [nlc], options)
    assert all(processed_nlc.lb == [-np.inf, -np.inf])
    assert all(processed_nlc.ub == [0, 0])
    assert options['nlconstr0'] == [2.1, 2.1]
