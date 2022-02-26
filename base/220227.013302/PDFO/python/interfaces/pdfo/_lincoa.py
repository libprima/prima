# -*- coding: utf-8 -*-
import warnings
from inspect import stack

import numpy as np


def lincoa(fun, x0, args=(), bounds=None, constraints=(), options=None):
    """LINCOA: LINearly Constrained Optimization Algorithm

    M. J. D. Powell did not publish any paper introducing LINCOA.

    Parameters
    ----------
    fun: callable
        The objective function, which accepts a vector `x` at input and returns a scalar.
    x0: ndarray, shape (n,)
        The initial guess. The size of `x0` should be consistent with the objective function.
    args: tuple, optional
        The extra-arguments to pass to the objective function. For example,

            ``lincoa(fun, x0, args, ...)``

        is equivalent to

            ``lincoa(lambda x: fun(x, *args), x0, ...)``

    bounds: ndarray of tuple with shape(n,2), or Bounds, optional
        Bound constraints of the problem. It can be one of the two cases below.
            1. An ndarray with shape(n,2). If the ndarray is 'bounds', then the bound constraint for x[i] is
                bounds[i, 0]<=x[i]<=bounds[i, 1]. Set bounds[i, 0] to -numpy.inf or None if there is no lower bound, and
                set bounds[i, 1] to numpy.inf or None if there is no upper bound.
            2. An instance of the `Bounds` class. Bounds(lb, ub) specifies a bound constraint lb<=x<=ub.
    constraints: LinearConstraint or a list of them, optional
        Linear constraints of the problem. It can be one of the two cases below.
            1. An instance of the `LinearConstraint` class.
                LinearConstraint(A, lb, ub) specifies a linear constraint lb<=A*x<=ub.
            2. A list of LinearConstraint.
    options: dict, optional
        The options passed to the solver. It is a structure that contains optionally:
            rhobeg: float, optional
                Initial value of the trust region radius, which should be a positive scalar. Typically,
                `options['rhobeg']` should be typically in the order of one tenth of the greatest expected change to a
                variable. By default, it is 1 if the problem is not scaled, 0.5 if the problem is scaled.
            rhoend: float, optional
                Final value of the trust region radius, which should be a positive scalar. `options['rhoend']` should
                indicate typically the accuracy required in the final values of the variables. Moreover,
                `options['rhoend']` should be no more than `options['rhobeg']` and is by default 1e-6.
            maxfev: int, optional
                Upper bound of the number of calls of the objective function `fun`. Its value must be not less than
                `options['npt']`+1. By default, it is 500*n.
            npt: int, optional
                Number of interpolation points of each model used in Powell's Fortran code. By default, it is 2*n+1.
            ftarget: float, optional
                Target value of the objective function. If a feasible iterate achieves an objective function value lower
                or equal to `options['ftarget']`, the algorithm stops immediately. By default, it is -numpy.inf.
            scale: bool, optional
                Flag indicating whether to scale the problem according to the bound constraints or not. By default, it
                is False. If the problem is to be scaled, then rhobeg and rhoend mentioned above will be used as the
                initial and final trust region radii for the scaled problem.
            quiet: bool, optional
                Flag of quietness of the interface. If it is set to True, the output message will not be printed. This
                flag does not interfere with the warning and error printing.
            classical: bool, optional
                Flag indicating whether to call the classical Powell code or not. By default, it is False.
            eliminate_lin_eq: bool, optional
                Flag indicating whether the linear equality constraints should be eliminated. By default, it is True.
            debug: bool, optional
                Debugging flag. By default, it is False.
            chkfunval: bool, optional
                Flag used when debugging. If both `options['debug']` and `options['chkfunval']` are True, an extra
                function evaluation would be performed to check whether the returned objective function value matches
                the returned x. By default, it is False.

    Returns
    -------
    res: OptimizeResult
        The results of the solver. Check ``pdfo.OptimizeResult`` for a description of the attributes.

    Notes
    -----
    The signature of this function is consistent with the `minimize` function available in ``scipy.optimize``, included
    in the SciPy package.

    See https://www.pdfo.net for more information.

    See also
    --------
    bobyqa : Bounded Optimization BY Quadratic Approximations
    cobyla : Constrained Optimization BY Linear Approximations
    newuoa : NEW Unconstrained Optimization Algorithm
    uobyqa : Unconstrained Optimization BY Quadratic Approximation
    pdfo : Powell's Derivative-Free Optimization solvers

    Examples
    --------
    1. The following code

    >>> from pdfo import *
    >>> import numpy as np
    >>> lin_con = LinearConstraint(2, None, 3)
    >>> options = {'maxfev': 50}
    >>> lincoa(np.cos, -1, constraints=lin_con, options=options)

    solves
        min  cos(x)
        s.t. 2*x <= 3
    starting from x0 = -1 with at most 50 function evaluations.

    2. The following code

    >>> from pdfo import *
    >>> obj = lambda x: x[0]**2 + x[1]**2
    >>> bounds = Bounds([0, 0.5], [2, 3])
    >>> lin_con = LinearConstraint([1, 1], 0, 1)
    >>> options = {'maxfev': 200}
    >>> lincoa(obj, [0, 1], bounds=bounds, constraints=lin_con, options=options)

    solves
        min  x^2 + y^2
        s.t. 0 <= x <= 2
             0.5 <= y <= 3
             0 <= x + y <= 1
    starting from [x0, y0] = [0, 1] with at most 200 function evaluations.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    try:
        from .gethuge import gethuge
    except ImportError:
        from ._dependencies import import_error_so

        # If gethuge cannot be imported, the execution should stop because the package is most likely not built.
        import_error_so('gethuge')

    from ._dependencies import prepdfo, _augmented_linear_constraint, postpdfo

    fun_name = stack()[0][3]  # name of the current function
    if len(stack()) >= 3:
        invoker = stack()[1][3].lower()
    else:
        invoker = ''

    # A cell that records all the warnings.
    # Why do we record the warning message in output['warnings'] instead of prob_info['warnings']? Because, if lincoa is
    # called by pdfo, then prob_info will not be passed to postpdfo, and hence the warning message will be lost. To the
    # contrary, output will be passed to postpdfo anyway.
    output = dict()
    output['warnings'] = []

    # Preprocess the inputs.
    fun_c, x0_c, bounds_c, constraints_c, options_c, _, prob_info = \
        prepdfo(fun, x0, args, bounds=bounds, constraints=constraints, options=options)

    # Check whether nonlinear constraints are passed to the function.
    if constraints_c['nonlinear'] is not None:
        warn_message = '{}: Nonlinear constraints are given as parameter; they will be ignored.'.format(fun_name)
        warnings.warn(warn_message, Warning)
        output['warnings'].append(warn_message)

    if invoker != 'pdfo' and prob_info['infeasible']:
        # The problem turned out infeasible during prepdfo.
        exitflag = -4
        nf = 1
        x = x0_c
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_x0']
        chist = np.array([constrviolation], dtype=np.float64)
        output['constr_modified'] = False
    elif invoker != 'pdfo' and prob_info['nofreex']:
        # x was fixed by the bound constraints during prepdfo.
        exitflag = 13
        nf = 1
        x = prob_info['fixedx_value']
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_fixedx']
        chist = np.array([constrviolation], dtype=np.float64)
        output['constr_modified'] = False
    elif invoker != 'pdfo' and prob_info['feasibility_problem']:
        # We could set fx=[], funcCount=0, and fhist=[] since no function evaluation occurred. But then we will have to
        # modify the validation of fx, funcCount, and fhist in postpdfo. To avoid such a modification, we set fx,
        # funcCount, and fhist as below and then revise them in postpdfo.
        nf = 1
        x = x0_c  # prepdfo has tried to set x0 to a feasible point (but may have failed)
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_x0']
        chist = np.array([constrviolation], dtype=np.float64)
        output['constr_value'] = np.asarray([], dtype=np.float64)
        if constrviolation < np.finfo(np.float64).eps:
            # Did prepdfo find a feasible point?
            exitflag = 14
        else:
            exitflag = 15
        output['constr_modified'] = False
    else:
        # The problem turns out 'normal' during prepdfo include all the constraints into one single linear constraint
        # (A_aug)'*x <= b_aug; note the TRANSPOSE due to the data structure of the Fortran code.
        n = x0_c.size
        a_aug, b_aug = _augmented_linear_constraint(n, bounds_c, constraints_c)
        a_aug = a_aug.T

        # Extract the options and parameters
        npt = options_c['npt']
        maxfev = options_c['maxfev']
        rhobeg = options_c['rhobeg']
        rhoend = options_c['rhoend']
        ftarget = options_c['ftarget']

        # The largest integer in the fortran functions; the factor 0.99 provides a buffer.
        max_int = np.floor(0.99 * gethuge('integer'))
        m = b_aug.size  # linear constraints: A_aug.T * x <= b_aug

        # The smallest nw, i.e., the nw with npt = n + 2. If it is larger than a threshold (system dependent), the
        # problem is too large to be executed on the system.
        min_nw = m * (2 + n) + (n + 2) * (2 * n + 6) + n * (9 + 3 * n) + max(m + 3 * n, 2 * m + n, 2 * n + 4)
        if min_nw >= max_int:
            executor = invoker.lower() if invoker == 'pdfo' else fun_name
            # nw would suffer from overflow in the Fortran code, exit immediately.
            raise SystemError('{}: problem too large for {}. Try other solvers.'.format(executor, fun_name))

        # The largest possible value for npt given that nw <= max_int.
        alpha = n + 7
        beta = 2 * m + m * (2 + n) + n * (9 + 3 * n) - max_int
        max_npt = max(n + 2, np.floor(0.5 * (-alpha + np.sqrt(alpha * alpha - 4 * beta))))
        if npt > max_npt:
            npt = max_npt
            w_message = \
                '{}: npt is so large that it is unable to allocate the workspace; it is set to {}'.format(fun_name, npt)
            warnings.warn(w_message, Warning)
            output['warnings'].append(w_message)
        if maxfev > max_int:
            maxfev = max_int
            w_message = \
                '{}: maxfev exceeds the upper limit of Fortran integer; it is set to {}'.format(fun_name, maxfev)
            warnings.warn(w_message, Warning)
            output['warnings'].append(w_message)

        # If x0 is not feasible, LINCOA will modify the constraints to make it feasible (which is a bit strange).
        # prepdfo has tried to find a feasible x0. Raise a warning is x0 is not 'feasible enough' so that the
        # constraints will be modified.
        if a_aug.size > 0 and any(np.dot(x0_c.T, a_aug) > b_aug + 1e-10 * max(1, np.max(b_aug))):
            output['constr_modified'] = True
            w_message = \
                '{}: preprocessing code did not find a feasible x0; problem is likely infeasible or SciPy is not ' \
                'installed on the machine; {} will modify the right-hand side of the constraints to make x0 ' \
                'feasible.'.format(fun_name, fun_name)
            warnings.warn(w_message, Warning)
            output['warnings'].append(w_message)
        else:
            output['constr_modified'] = False

        # Call the Fortran code.
        try:
            if options_c['classical']:
                from . import flincoa_classical as flincoa
            else:
                from . import flincoa
        except ImportError:
            from ._dependencies import import_error_so
            import_error_so()

        # m should be precised not to raise any error if there is no linear constraints.
        x, fx, exitflag, fhist, chist, constrviolation = \
            flincoa.mlincoa(npt, m, a_aug, b_aug, x0_c, rhobeg, rhoend, 0, maxfev, ftarget, fun_c)
        nf = int(flincoa.flincoa.nf)

    # Postprocess the result.
    return postpdfo(x, fx, exitflag, output, fun_name, nf, fhist, options_c, prob_info, constrviolation, chist)
