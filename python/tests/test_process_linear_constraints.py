import numpy as np
from prima import LinearConstraint, process_multiple_linear_constraints, separate_LC_into_eq_and_ineq


def test_multiple_linear_constraints():
    constraints = [LinearConstraint(A=np.array([[1, 2], [3, 4]]), lb=[5, 6], ub=[7, 8]),
                   LinearConstraint(A=np.array([[9, 10], [11, 12]]), lb=[13, 14], ub=[15, 16])]
    processed_constraint = process_multiple_linear_constraints(constraints)
    assert (processed_constraint.A == np.array([[1, 2], [3, 4], [9, 10], [11, 12]])).all()
    assert all(processed_constraint.lb == [5, 6, 13, 14])
    assert all(processed_constraint.ub == [7, 8, 15, 16])


def test_separate_LC_into_eq_and_ineq():
    linear_constraint = LinearConstraint(A=np.array([[1, 2], [3, 4]]), lb=[5, 6], ub=[5, 8])
    A_eq, b_eq, A_ineq, b_ineq = separate_LC_into_eq_and_ineq(linear_constraint)
    assert (A_eq == np.array([[1, 2]])).all()
    assert all(b_eq == [5])
    assert (A_ineq == np.array([[3, 4], [-3, -4]])).all()
    assert all(b_ineq == [8, -6])
