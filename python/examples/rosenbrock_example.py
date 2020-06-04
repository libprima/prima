#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Illustration of how to use pdfo.

Authors
-------
Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
Department of Applied Mathematics,
The Hong Kong Polytechnic University.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
"""
from pdfo import pdfo, Bounds, LinearConstraint, NonlinearConstraint
# If SciPy (version 1.1 or above) is installed, then Bounds, LinearConstraint,
# and NonlinearConstraint can alternatively be imported from scipy.optimize.
import numpy as np


def chrosen(x):  # the subroutine defining the objective function
    """Chained Rosenbrock function."""
    return sum((1 - x[:-1]) ** 2 + 4 * (x[1:] - x[:-1] ** 2) ** 2)


def nlc_ineq(x):  # the subroutine defining the nonlinear inequality constraints
    """Example of nonlinear inequality constraint function."""
    return -x[:-1] ** 2 + x[1:]


def nlc_eq(x):  # the subroutine defining the nonlinear equality constraints
    """Example of nonlinear equality constraint function (unit sphere)."""
    return sum(x ** 2) - 1


if __name__ == '__main__':
    print('\nMinimize the chained Rosenbrock function with three variables subject to various constraints:\n')
    np.set_printoptions(precision=4, threshold=7, edgeitems=3)  # printing NumPy arrays
    x0 = [0, 0, 0]  # starting point

    print('\n1. Nonlinear constraints --- ||x||_2^2 = 1, x(i)^2 >= x(i+1) >= 0.5*x(i) >= 0 for i = 1, 2:\n')
    # bound constraints lb <= x <= ub
    lb = [0, 0, 0]
    ub = [np.inf, np.inf, np.inf]  # ub = [None, None, None] works equally well
    bounds = Bounds(lb, ub)  # bounds = [(lb[0], ub[0]), (lb[1], ub[1]), (lb[2], ub[2])] works equally well
    # linear inequality constraints lin_lb <= A*x <= lin_ub
    A = [[0.5, -1, 0], [0, 0.5, -1]]
    lin_lb = [-np.inf, -np.inf]
    lin_ub = [0, 0]
    lin_con = LinearConstraint(A, lin_lb, lin_ub)
    # nonlinear inequality constraints: nonlin_lb <= nlc_ineq(x) <= nonlin_ub
    nonlin_lb = [-np.inf, -np.inf]
    nonlin_ub = [0, 0]
    nonlin_con_ineq = NonlinearConstraint(nlc_ineq, nonlin_lb, nonlin_ub)
    nonlin_con_eq = NonlinearConstraint(nlc_eq, 0, 0)  # this is an equality constraint
    constraints = [lin_con, nonlin_con_ineq, nonlin_con_eq]  # all constraints
    res = pdfo(chrosen, x0, bounds=bounds, constraints=constraints)
    print(res)
    # The following representation of the constraints works also well
    # con = lambda x: np.concatenate((np.dot(A, x), nlc_ineq(x), [nlc_eq(x)]))
    # con_lb = np.concatenate((lin_lb, nonlin_lb, [0]))
    # con_ub = np.concatenate((lin_ub, nonlin_ub, [0]))
    # constraints = NonlinearConstraint(con, con_lb, con_ub)

    print('\n2. Linear constraints --- sum(x) = 1, x(i+1) <= x(i) <= 1 for i = 1, 2:\n')
    bounds = Bounds([-np.inf, -np.inf, -np.inf], [1, 1, 1])
    A = [[-1, 1, 0], [0, -1, 1], [1, 1, 1]]
    lin_con = LinearConstraint(A, [-np.inf, -np.inf, 1], [0, 0, 1])
    res = pdfo(chrosen, x0, bounds=bounds, constraints=lin_con)
    print(res)

    print('\n3. Bound constraints --- -0.5 <= x(1) <= 0.5, 0 <= x(2) <= 0.25:\n')
    bounds = Bounds([-0.5, 0, -np.inf], [0.5, 0.25, np.inf])
    options = {'rhobeg': 0.1}
    res = pdfo(chrosen, x0, bounds=bounds, options=options)
    print(res)

    print('\n4. No constraints:\n')
    res = pdfo(chrosen, x0)
    print(res)
