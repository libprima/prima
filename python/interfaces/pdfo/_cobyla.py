# -*- coding: utf-8 -*-
import warnings
from inspect import stack

import numpy as np


def cobyla(fun, x0, args=(), bounds=None, constraints=(), options=None):
    """COBYLA: Constrained Optimization BY Linear Approximations

    The algorithm is described in [M. J. D. Powell, A direct search optimization method that models the objective and
    constraint functions by linear interpolation, In Advances in Optimization and Numerical Analysis, eds. S. Gomez and
    J. P. Hennart, pages 51--67, Springer Verlag, Dordrecht, Netherlands, 1994].

    Parameters
    ----------
    fun: callable
        The objective function, which accepts a vector `x` at input and returns a scalar.
    x0: ndarray, shape (n,)
        The initial guess. The size of `x0` should be consistent with the objective function.
    args: tuple, optional
        The extra-arguments to pass to the objective function. For example,

            ``cobyla(fun, x0, args, ...)``

        is equivalent to

            ``cobyla(lambda x: fun(x, args), x0, ...)``

    bounds: either ndarray of tuple with shape(n,2), or Bounds, optional
        Bound constraints of the problem. The bounds can be specified in two different ways:
            1. Instance of `Bounds` class.
            2. Sequence of (lb, ub) pairs for each element in `x`. To specify that `x[i]` is unbounded below, set
               `bounds[i, 0]` to -np.inf; set `bounds[i, 1]` to np.inf if `x[i]` is unbounded above.
    constraints: dict, LinearConstraint, NonlinearConstraint or list of them, optional
        Constraints of the problem. It can be one of the three cases below.
            1. A dictionary with fields:
                type: str
                    Constraint type: 'eq' for equality constraints and 'ineq' for inequality constraints.
                fun: callable
                    The constraint function.
            2. Instances of LinearConstraint or NonlinearConstraint.
            3. A list, each of whose elements can be a dictionary described in 1, an instance of LinearConstraint or an
               instance of NonlinearConstraint.
    options: dict, optional
        The options passed to the solver. It is a structure that contains optionally:
            rhobeg: float, optional
                Initial value of the trust region radius, which should be a positive scalar. `options['rhobeg']` should
                be typically set roughly to one tenth of the greatest expected change to a variable. By default, it is
                1.
            rhoend: float, optional
                Final value of the trust region radius, which should be a positive scalar. `options['rhoend']` should
                indicate typically the accuracy required in the final values of the variables. Moreover,
                `options['rhoend']` should be no more than `options['rhobeg']` and is by default 1e-6.
            maxfev: int, optional
                Upper bound of the number of calls of the objective function `fun`. Its value must be not less than
                n+2. By default, it is 500*n.
            ftarget: float, optional
                Target value of the objective function. If a feasible iterate achieves an objective function value lower
                or equal to `options['ftarget']`, the algorithm stops immediately. By default, it is -np.inf.
            scale: bool, optional
                Flag indicating whether to scale the problem. If it is True, the variables will be scaled according to
                the bounds constraints if any. By default, it is False.
            quiet: bool, optional
                Flag of quietness of the interface. If it is set to True, the output message will not be printed. This
                flag does not interfere with the warning and error printing.
            classical: bool, optional
                Flag indicating whether to call the classical Powell code or not. By default, it is False.
            debug: bool, optional
                Debugging flag. By default, it is False.
            chkfunval: bool, optional
                Flag used when debugging. If both `options['debug']` and `options['chkfunval']` are True, an extra
                function evaluation would be performed to check whether the returned objective function value is
                consistent with the returned x. By default, it is False.

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
    lincoa : LINearly Constrained Optimization Algorithm
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
    >>> cobyla(np.cos, -1, constraints=lin_con, options=options)

    solves
        min  cos(x)
        s.t. 2*x <= 3
    starting from x0 = -1 with at most 50 function evaluations.

    2. The following code

    >>> from pdfo import *
    >>> obj = lambda x: x[0]**2 + x[1]**2
    >>> con = lambda x: x[0]**2 - x[1]
    >>> bounds = Bounds([0, 0.5], [2, 3])
    >>> lin_con = LinearConstraint([1, 1], 0, 1)
    >>> nonlin_con = NonlinearConstraint(con, None, 0)
    >>> options = {'maxfev': 200}
    >>> cobyla(obj, [0, 1], bounds=bounds, constraints=[lin_con, nonlin_con], options=options)

    solves
        min  x^2 + y^2
        s.t. 0 <= x <= 2
             0.5 <= y <= 3
             0 <= x + y <= 1
             x^2 - y <= 0
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
    # Why do we record the warning message in output['warnings'] instead of prob_info['warnings']? Because, if cobyla is
    # called by pdfo, then prob_info will not be passed to postpdfo, and hence the warning message will be lost. To the
    # contrary, output will be passed to postpdfo anyway.
    output = dict()
    output['warnings'] = []

    # Preprocess the inputs.
    fun_c, x0_c, bounds_c, constraints_c, options_c, _, prob_info = \
        prepdfo(fun, x0, args, bounds=bounds, constraints=constraints, options=options)

    if prob_info['infeasible']:
        # The problem turned out infeasible during prepdfo.
        exitflag = -4
        nf = 1
        x = x0_c
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_x0']
        chist = np.array([constrviolation], dtype=np.float64)
        output['nlc'] = prob_info['nlc_x0']
    elif prob_info['nofreex']:
        # x was fixed by the bound constraints during prepdfo
        exitflag = 13
        nf = 1
        x = prob_info['fixedx_value']
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_fixedx']
        chist = np.array([constrviolation], dtype=np.float64)
        output['nlc'] = prob_info['nlc_fixedx']
    else:
        # The problem turns out 'normal' during prepdfo include all the constraints into one single nonlinear
        # constraint.
        n = x0_c.size
        a_aug, b_aug = _augmented_linear_constraint(n, bounds_c, constraints_c)

        def ctr(x_aug):
            c = np.array([], dtype=np.float64)

            if b_aug.size > 0:
                cx = np.dot(a_aug, x_aug) - b_aug
                c = np.concatenate((c, cx))

            if constraints_c['nonlinear'] is not None:
                cx = constraints_c['nonlinear']['fun'](x_aug)
                c = np.concatenate((c, cx))

            return -c

        # The constraint function should be evaluated only one time for each array, not one time for each element of the
        # array. Since the Python list is a mutable structure, it is used to create a memory, by never changing the
        # reference to it in the constraint function.
        x_aug_mem = []

        def ctr_elmt(x_aug, i):
            if i == 0:
                con_x = ctr(x_aug)
                del x_aug_mem[:]  # the old memory is removed without changing the reference to the list
                x_aug_mem.extend(con_x)
                return con_x[i]
            else:
                return x_aug_mem[i]

        conval_x0 = ctr(x0_c)
        m = conval_x0.size

        # Extract the options and parameters.
        maxfev = options_c['maxfev']
        rhobeg = options_c['rhobeg']
        rhoend = options_c['rhoend']
        ftarget = options_c['ftarget']

        # The largest integer in the fortran functions; the factor 0.99 provides a buffer.
        max_int = np.floor(0.99 * gethuge('integer'))

        # The smallest nw, i.e., the nw with npt = n + 2. If it is larger than a threshold (system dependent), the
        # problem is too large to be executed on the system.
        min_nw = n * (3 * n + 2 * m + 11) + 4 * m + 6
        if min_nw >= max_int:
            executor = invoker.lower() if invoker == 'pdfo' else fun_name
            # nw would suffer from overflow in the Fortran code, exit immediately.
            raise SystemError('{}: problem too large for {}. Try other '
                              'solvers.'.format(executor, fun_name))
        if maxfev > max_int:
            maxfev = max_int
            w_message = '{}: maxfev exceeds the upper limit of Fortran integer; it is set to ' \
                        '{}'.format(fun_name, maxfev)
            warnings.warn(w_message, Warning)
            output['warnings'].append(w_message)

        # Call the Fortran code.
        try:
            if options_c['classical']:
                from . import fcobyla_classical as fcobyla
            else:
                from . import fcobyla
        except ImportError:
            from ._dependencies import import_error_so
            import_error_so()

        # m should be precised not to raise any error if there is no linear constraints.
        x, fx, exitflag, fhist, chist, constrviolation, conval = \
            fcobyla.mcobyla(x0_c, rhobeg, rhoend, 0, maxfev, ftarget, conval_x0, fun_c,
                            lambda i_c, x_c: ctr_elmt(x_c, i_c - 1))
        nf = int(fcobyla.fcobyla.nf)

        if m > 0:
            output['nlc'] = -conval[b_aug.size:]

    # Postprocess the result.
    return postpdfo(x, fx, exitflag, output, fun_name, nf, fhist, options_c, prob_info, constrviolation, chist)
