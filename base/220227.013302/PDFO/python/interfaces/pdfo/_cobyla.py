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

            ``cobyla(lambda x: fun(x, *args), x0, ...)``

    bounds: ndarray of tuple with shape(n,2), or Bounds, optional
        Bound constraints of the problem. It can be one of the two cases below. 
            1. An ndarray with shape(n,2). If the ndarray is 'bounds', then the bound constraint for x[i] is 
                bounds[i, 0]<=x[i]<=bounds[i, 1]. Set bounds[i, 0] to -numpy.inf or None if there is no lower bound, and 
                set bounds[i, 1] to numpy.inf or None if there is no upper bound. 
            2. An instance of the `Bounds` class. Bounds(lb, ub) specifies a bound constraint lb<=x<=ub.
    constraints: dict, LinearConstraint, NonlinearConstraint, or a list of them, optional
        Constraints of the problem. It can be one of the three cases below.
            1. A dictionary with fields:
                type: str
                    Constraint type: 'eq' for equality constraints and 'ineq' for inequality constraints.
                fun: callable
                    The constraint function.
                When type='eq', such a dictionary specifies an equality constraint fun(x)=0; 
                when type='ineq', it specifies an inequality constraint fun(x)>=0.
            2. An instance of the `LinearConstraint` class or the `NonlinearConstraint` class.
                LinearConstraint(A, lb, ub) specifies a linear constraint lb<=A*x<=ub;
                NonLinearConstraint(fun, lb, ub) specifies a nonlinear constraint lb<=fun(x)<=ub.
    options: dict, optional
        The options passed to the solver. It is a structure that contains optionally:
            rhobeg: float, optional
                Initial value of the trust region radius, which should be a positive scalar. Typically, `options['rhobeg']` 
                should be in the order of one tenth of the greatest expected change to a variable. By default, it is 1 if 
                the problem is not scaled, 0.5 if the problem is scaled.
            rhoend: float, optional
                Final value of the trust region radius, which should be a positive scalar. `options['rhoend']` should
                indicate typically the accuracy required in the final values of the variables. Moreover,
                `options['rhoend']` should be no more than `options['rhobeg']` and is by default 1e-6.
            maxfev: int, optional
                Upper bound of the number of calls of the objective function `fun`. Its value must be not less than
                n+2. By default, it is 500*n.
            ftarget: float, optional
                Target value of the objective function. If a feasible iterate achieves an objective function value lower
                or equal to `options['ftarget']`, the algorithm stops immediately. By default, it is -numpy.inf.
            scale: bool, optional
                Flag indicating whether to scale the problem acording to the bound constraints. By default, it is False. 
                If the problem is to be scaled, then rhobeg and rhoend mentioned above will be used as the initial and 
                final trust-region radii for the scaled problem.
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
                function/constriant evaluation would be performed to check whether the returned values of the objective 
                function and the constraint match the returned x. By default, it is False.

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

    if invoker != 'pdfo' and prob_info['infeasible']:
        # The problem turned out infeasible during prepdfo.
        exitflag = -4
        nf = 1
        x = x0_c
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_x0']
        chist = np.array([constrviolation], dtype=np.float64)
        output['constr_value'] = prob_info['nlc_x0']
    elif invoker != 'pdfo' and prob_info['nofreex']:
        # x was fixed by the bound constraints during prepdfo
        exitflag = 13
        nf = 1
        x = prob_info['fixedx_value']
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_fixedx']
        chist = np.array([constrviolation], dtype=np.float64)
        output['constr_value'] = prob_info['nlc_fixedx']
    elif invoker != 'pdfo' and prob_info['feasibility_problem'] and \
            prob_info['refined_type'] != 'nonlinearly-constrained':
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
    else:
        # The problem turns out 'normal' during prepdfo include all the constraints into one single nonlinear
        # constraint.
        n = x0_c.size
        a_aug, b_aug = _augmented_linear_constraint(n, bounds_c, constraints_c)

        # The constraint function received by COBYLA can return an array: in fact, the Fortran code interpret this
        # function as a subroutine from v1.0.
        def ctr(x_aug):
            c = np.array([], dtype=np.float64)

            if b_aug.size > 0:
                cx = np.dot(a_aug, x_aug) - b_aug
                c = np.concatenate((c, cx))

            if constraints_c['nonlinear'] is not None:
                cx = constraints_c['nonlinear']['fun'](x_aug)
                c = np.concatenate((c, cx))

            return -c

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
            fcobyla.mcobyla(x0_c, rhobeg, rhoend, 0, maxfev, ftarget, conval_x0, fun_c, lambda m, x: ctr(x))
        nf = int(fcobyla.fcobyla.nf)

        if m > 0:
            output['constr_value'] = -conval[b_aug.size:]

    # Postprocess the result.
    return postpdfo(x, fx, exitflag, output, fun_name, nf, fhist, options_c, prob_info, constrviolation, chist)
