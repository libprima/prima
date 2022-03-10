# -*- coding: utf-8 -*-
import numpy as np


def pdfo(fun, x0, args=(), method=None, bounds=None, constraints=(), options=None):
    """PDFO: Powell's Derivative-Free Optimization solvers

    PDFO provides an interface to call Powell's derivatives-free optimization solvers: UOBYQA, NEWUOA, BOBYQA, LINCOA,
    COBYLA.

    UOBYQA: Unconstrained Optimization BY Quadratic Approximation
    The algorithm is described in [M. J. D. Powell, UOBYQA: unconstrained optimization by quadratic approximation, Math.
    Program., 92(B):555--582, 2002].

    NEWUOA: NEW Unconstrained Optimization Algorithm
    The algorithm is described in [M. J. D. Powell, The NEWUOA software for unconstrained optimization without
    derivatives, In Large-Scale Nonlinear Optimization, eds. G. Di Pillo and M. Roma, pages 255--297, Springer, New
    York, US, 2006].

    BOBYQA: Bounded Optimization BY Quadratic Approximation
    The algorithm is described in [M. J. D. Powell. The BOBYQA algorithm for bound constrained optimization without
    derivatives, Technical Report DAMTP 2009/NA06, Department of Applied Mathematics and Theoretical Physics, Cambridge
    University, Cambridge, UK, 2009].

    LINCOA: LINearly Constrained Optimization Algorithm
    M. J. D. Powell did not publish any paper introducing LINCOA.

    COBYLA: Constrained Optimization BY Linear Approximations
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

            ``pdfo(fun, x0, args, ...)``

        is equivalent to

            ``pdfo(lambda x: fun(x, *args), x0, ...)``

    method: str, optional
        The name of the Powell method that will be used.
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
            3. A list, each of whose elements is a dictionary described in 1, a LinearConstraint, or a
                NonlinearConstraint.
    options: dict, optional
        The options passed to the solver. It is a structure that contains optionally:
            rhobeg: float, optional
                Initial value of the trust region radius, which should be a positive scalar. Typically,
                `options['rhobeg']` should be in the order of one tenth of the greatest expected change to a variable.
                By default, it is 1 if the problem is not scaled (but min(1, min(ub-lb)/4) if the solver is BOBYQA), 0.5
                if the problem is scaled.
            rhoend: float, optional
                Final value of the trust region radius, which should be a positive scalar. `options['rhoend']` should
                indicate the accuracy required in the final values of the variables. Moreover, `options['rhoend']`
                should be no more than `options['rhobeg']` and is by default 1e-6.
            maxfev: int, optional
                Upper bound of the number of calls of the objective function `fun`. Its value must be not less than
                `options['npt']`+1. By default, it is 500*n.
            npt: int, optional
                Number of interpolation points of each model used in Powell's Fortran code. It is used only if the
                solver is NEWUOA, BOBYQA, or LINCOA.
            ftarget: float, optional
                Target value of the objective function. If a feasible iterate achieves an objective function value lower
                or equal to `options['ftarget']`, the algorithm stops immediately. By default, it is -numpy.inf.
            scale: bool, optional
                Flag indicating whether to scale the problem according to the bound constraints. By default, it is
                False. If the problem is to be scaled, then rhobeg and rhoend mentioned above will be used as the
                initial and final trust region radii for the scaled problem.
            honour_x0: bool, optional
                Flag indicating whether to respect the user-defined x0. By default, it is False. It is used only if the 
                solver is BOBYQA.
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
                function/constraint evaluation would be performed to check whether the returned values of objective
                function and constraint match the returned x. By default, it is False.

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
    uobyqa : Unconstrained Optimization BY Quadratic Approximation
    newuoa : NEW Unconstrained Optimization Algorithm
    bobyqa : Bounded Optimization BY Quadratic Approximations
    lincoa : LINearly Constrained Optimization Algorithm
    cobyla : Constrained Optimization BY Linear Approximations

    Examples
    --------
    1. The following code

    >>> from pdfo import *
    >>> import numpy as np
    >>> lin_con = LinearConstraint(2, None, 3)
    >>> options = {'maxfev': 50}
    >>> pdfo(np.cos, -1, constraints=lin_con, options=options)

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
    >>> pdfo(obj, [0, 1], bounds=bounds, constraints=[lin_con, nonlin_con], options=options)

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
    from ._dependencies import prepdfo, postpdfo

    # A cell that records all the warnings.
    output = dict()
    output['warnings'] = []

    # Preprocess the inputs.
    fun_c, x0_c, bounds_c, constraints_c, options_c, method, prob_info = \
        prepdfo(fun, x0, args, method, bounds, constraints, options)

    if prob_info['infeasible']:
        # The problem turned out infeasible during prepdfo.
        exitflag = -4
        nf = 1
        x = x0_c
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_x0']
        chist = np.array([constrviolation], dtype=np.float64)
        output['constr_value'] = prob_info['nlc_x0']
        output['constr_modified'] = False
    elif prob_info['nofreex']:
        # x was fixed by the bound constraints during prepdfo.
        exitflag = 13
        nf = 1
        x = prob_info['fixedx_value']
        fx = fun_c(x)
        fhist = np.array([fx], dtype=np.float64)
        constrviolation = prob_info['constrv_fixedx']
        chist = np.array([constrviolation], dtype=np.float64)
        output['constr_value'] = prob_info['nlc_fixedx']
        output['constr_modified'] = False
    elif prob_info['feasibility_problem'] and prob_info['refined_type'] != 'nonlinearly-constrained':
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
        # The problem turns out 'normal' during prepdfo.
        lower_method = method.lower()
        try:
            if lower_method == 'uobyqa':
                from . import uobyqa
                opti_res = uobyqa(fun_c, x0_c, options=options_c)
            elif lower_method == 'newuoa':
                from . import newuoa
                opti_res = newuoa(fun_c, x0_c, options=options_c)
            elif lower_method == 'bobyqa':
                from . import bobyqa
                opti_res = bobyqa(fun_c, x0_c, bounds=bounds_c, options=options_c)
            elif lower_method == 'lincoa':
                from . import lincoa
                opti_res = lincoa(fun_c, x0_c, bounds=bounds_c, constraints=constraints_c, options=options_c)
            elif lower_method == 'cobyla':
                from . import cobyla
                opti_res = cobyla(fun_c, x0_c, bounds=bounds_c, constraints=constraints_c, options=options_c)
        except ImportError:
            from ._dependencies import import_error_so
            import_error_so(lower_method)

        # Extract the output from the solvers. The output is extended with the possible outputs returned by some
        # specific solvers (like the nonlinear constraint nlc for COBYLA).
        x = opti_res.x
        fx = opti_res.fun
        exitflag = opti_res.status
        nf = opti_res.nfev
        fhist = opti_res.fhist
        try:
            constrviolation = opti_res.constrviolation
        except AttributeError:
            constrviolation = 0
        try:
            chist = opti_res.chist
        except AttributeError:
            chist = None
        try:
            output['constr_value'] = opti_res.constr_value
        except AttributeError:
            pass
        try:
            output['constr_modified'] = opti_res.constr_modified
        except AttributeError:
            pass

        # The warnings that have been raised in the solvers and treated during their own calls to postpdfo should be
        # transfer to the call to postpdfo of pdfo to appear to the output of pdfo.
        output['warnings'].extend(opti_res.warnings)

    # Postprocess the result.
    return postpdfo(x, fx, exitflag, output, method, nf, fhist, options_c, prob_info, constrviolation, chist)
