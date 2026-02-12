# Bounds may appear unused in this file but we need to import it to make it available to the user
from scipy.optimize import NonlinearConstraint, LinearConstraint, Bounds, OptimizeResult
from ._nonlinear_constraints import process_nl_constraints
from ._linear_constraints import (
    combine_multiple_linear_constraints,
    separate_LC_into_eq_and_ineq,
)
from ._bounds import process_bounds
from enum import Enum
import numpy as np
from ._common import _project
from ._common import get_arrays_tol
from .infos import FIXED_SUCCESS
from warnings import warn

# TODO: Set __version__ without going to the bindings

class PRIMAMessage(Enum):
    # See prima_message_t in prima.h
    NONE = 0
    EXIT = 1
    RHO = 2
    FEVL = 3


class ConstraintType(Enum):
    LINEAR_OBJECT = 5
    NONLINEAR_OBJECT = 10
    LINEAR_DICT = 15
    NONLINEAR_DICT = 20


def get_constraint_type(constraint):
    if isinstance(constraint, dict) and ("A" in constraint) and ("lb" in constraint) and ("ub" in constraint):
        return ConstraintType.LINEAR_DICT
    elif isinstance(constraint, dict) and ("fun" in constraint) and ("lb" in constraint) and ("ub" in constraint):
        return ConstraintType.NONLINEAR_DICT
    elif hasattr(constraint, "A") and hasattr(constraint, "lb") and hasattr(constraint, "ub"):
        return ConstraintType.LINEAR_OBJECT
    elif hasattr(constraint, "fun") and hasattr(constraint, "lb") and hasattr(constraint, "ub"):
        return ConstraintType.NONLINEAR_OBJECT
    else:
        raise ValueError("Constraint type not recognized")


def process_constraints(constraints):
    # First throw it back if it's an empty tuple
    if not constraints:
        return None, None
    # Next figure out if it's a list of constraints or a single constraint
    # If it's a single constraint, make it a list, and then the remaining logic
    # doesn't have to change
    if not isinstance(constraints, list):
        constraints = [constraints]

    # Separate out the linear and nonlinear constraints
    linear_constraints = []
    nonlinear_constraints = []
    for constraint in constraints:
        constraint_type = get_constraint_type(constraint)
        if constraint_type is ConstraintType.LINEAR_OBJECT:
            linear_constraints.append(constraint)
        elif constraint_type is ConstraintType.NONLINEAR_OBJECT:
            nonlinear_constraints.append(constraint)
        elif constraint_type == ConstraintType.LINEAR_DICT:
            linear_constraints.append(LinearConstraint(constraint["A"], constraint["lb"], constraint["ub"]))
        elif constraint_type == ConstraintType.NONLINEAR_DICT:
            nonlinear_constraints.append(NonlinearConstraint(constraint["fun"], constraint["lb"], constraint["ub"]))
        else:
            raise ValueError("Constraint type not recognized")

    if len(nonlinear_constraints) > 0:
        nonlinear_constraint_function = process_nl_constraints(nonlinear_constraints)
    else:
        nonlinear_constraint_function = None

    # Determine if we have multiple linear constraints, just 1, or none, and process accordingly
    if len(linear_constraints) > 1:
        linear_constraint = combine_multiple_linear_constraints(linear_constraints)
    elif len(linear_constraints) == 1:
        linear_constraint = linear_constraints[0]
    else:
        linear_constraint = None

    return linear_constraint, nonlinear_constraint_function


def minimize(fun, x0, args=(), method=None, bounds=None, constraints=(), callback=None, options=None):
    r'''Powell Reference Implementation for Modernization and Amelioration

    PRIMA is an interface to call Powell's derivatives-free optimization solvers:
    UOBYQA, NEWUOA, BOBYQA, LINCOA, and COBYLA. They are designed to minimize a
    scalar function of several variables subject to (possibly) bound
    constraints, linear constraints, and nonlinear constraints.

    PRIMA is dedicated to the late Professor M. J. D. Powell, FRS (1936-2015)

    It is written by Zaikun Zhang.
    
    Python bindings and Python implementation of COBYLA contributed by
    Nickolai Belakovski.

    Parameters
    ----------
    fun : callable
        Objective function to be minimized.

            ``fun(x, *args) -> float``

        where ``x`` is an array with shape (n,) and `args` is a tuple.
    
    x0 : array_like, shape (n,)
        Initial guess.

    args : tuple, optional
        Extra arguments of the objective function. For example,

            ``pdfo(fun, x0, args, ...)``

        is equivalent to

            ``pdfo(lambda x: fun(x, *args), x0, ...)``

    method : {'uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla'}, optional
        Name of the Powell method that will be used. By default, 'newuoa'
        is selected if the problem is unconstrained, 'bobyqa' is selected
        if the problem is bound-constrained, 'lincoa' is selected if the
        problem is linearly constrained, and 'cobyla' is selected if the
        problem is nonlinearly constrained.

    bounds : {`scipy.optimize.Bounds`, array_like, shape (n, 2)}, optional
        Bound constraints of the problem. It can be one of the cases below.

        #. An instance of `scipy.optimize.Bounds`.
        #. An array with shape (n, 2). The bound constraints for ``x[i]`` are
           ``bounds[i, 0] <= x[i] <= bounds[i, 1]``. Set ``bounds[i, 0]`` to
           :math:`-\infty` if there is no lower bound, and set ``bounds[i, 1]``
           to :math:`\infty` if there is no upper bound.

    constraints : {dict, `scipy.optimize.LinearConstraint`, `scipy.optimize.NonlinearConstraint`, list, tuple}, optional
        Constraints of the problem. It can be one of the cases below.

        #. A dictionary with fields:

            type : {'eq', 'ineq'}
                Whether the constraint is ``fun(x) = 0`` or ``fun(x) >= 0``.
            fun : callable
                Constraint function.

        #. An instance of `scipy.optimize.LinearConstraint`.
        #. An instance of `scipy.optimize.NonlinearConstraint`.
        #. A list or tuple, each of whose elements are described in 1, 2, and 3.

    options : dict, optional
        The options passed to the solver. Accepted options are:

            rhobeg : float, optional
                Reasonable initial changes to the variables.
            tol : float, optional
                Final accuracy in the optimization (not precisely guaranteed).
                This is a lower bound on the size of the trust region.
            iprint : int, optional
                Controls the frequency of output:
                    0. (default) There will be no printing
                    1. A message will be printed to the screen at the end of iteration, showing
                    the best vector of variables found and its objective function value
                    2. in addition to 1, each new value of RHO is printed to the screen,
                    with the best vector of variables so far and its objective function
                    value.
                    3. in addition to 2, each function evaluation with its variables will
                    be printed to the screen.
            maxfev : int, optional
                Maximum number of function evaluations.
            ctol : float, optional
                Tolerance (absolute) for constraint violations
            ftarget : float, optional
                Stop if the objective function is less than `ftarget`.

    Returns
    -------
    res : `scipy.optimize.OptimizeResult`
        Result of the optimization procedure, with the following fields:

            message : str
                Description of the exit status specified in the ``status``
                field (i.e., the cause of the termination of the solver).
            success : bool
                Whether the optimization procedure terminated successfully.
            status : int
                Termination status of the optimization procedure.
            fun : float
                Objective function value at the solution point.
            x : `numpy.ndarray`, shape (n,)
                Solution point.
            maxcv: float
                Maximum constraint violation at the solution point.
            nlconstr : `numpy.ndarray`
                The values of the nonlinear constraints at the solution point.
            nfev : int
                Number of function evaluations.
            method : str
                Name of the Powell method used.

        A description of the termination statuses is given below.

            .. list-table::
                :widths: 25 75
                :header-rows: 1

                * - Exit status
                - Description
                * - 0
                - The lower bound on the trust-region radius is reached.
                * - 1
                - The target value of the objective function is reached.
                * - 2
                - A trust-region step has failed to reduce the quadratic model.
                * - 3
                - The maximum number of function evaluations is reached.
                * - 6
                - There is not enough space between some lower and upper bounds, namely
                    one of the difference XU(I)-XL(I) is less than 2*RHOBEG
                * - 7
                - Rounding errors are becoming damaging.
                * - 8
                - One of the linear constraints has a zero gradient
                * - 13
                - All variables are fixed by the bounds.
                * - 20
                - The maximum number of iterations (trust region steps) is reached.
                * - 30
                - The iteration has been terminated by the callback function.
                * - -1
                - NaN is encountered in the solution point.
                * - -2
                - NaN or +Inf is encountered in the objective/constraint function value.
                * - -3
                - NaN is encountered in the model parameter.
                * - 100
                - Invalid input from the user.
                * - 101
                - The Fortran code reached an assertion.
                * - 102
                - The Fortran code reached a validation faillure (similar to assertion).
                * - 103
                - Memory allocation failure.

    References
    ----------
    Z. Zhang, PRIMA: Reference Implementation for Powell's Methods with Modernization and Amelioration,
        available at https://www.libprima.net, [DOI: 10.5281/zenodo.8052654](https://doi.org/10.5281/zenodo.8052654), 2023

    Examples
    ----------
    The following example shows how to solve a simple optimization problem using
    `prima`. In practice, the  problem considered below should be solved with a
    derivative-based method as it is a smooth problem for which the derivatives
    are known. We solve it here using `pdfo` only as an illustration.

    We consider the 2-dimensional problem

    .. math::

        \min_{x, y \in \R} \quad x^2 + y^2 \quad \text{s.t.} \quad \left\{
        \begin{array}{l}
            0 \le x \le 2,\\
            1 / 2 \le y \le 3,\\
            0 \le x + y \le 1,\\
            x^2 - y \le 0.
        \end{array} \right.

    We solve this problem using `prima` starting from the initial guess
    :math:`(x_0, y_0) = (0, 1)` with at most 200 function evaluations.

    .. testsetup::

        import numpy as np
        np.set_printoptions(precision=1, suppress=True)


    >>> import numpy as np
    >>> np.set_printoptions(precision=1, suppress=True)
    >>> from prima import minimize, Bounds, LinearConstraint, NonlinearConstraint
    >>>
    >>> # Build the constraints.
    >>> bounds = Bounds([0, 0.5], [2, 3])
    >>> linear_constraints = LinearConstraint([1, 1], 0, 1)
    >>> nonlinear_constraints = NonlinearConstraint(lambda x: x[0]**2 - x[1], -np.inf, 0)
    >>> constraints = [linear_constraints, nonlinear_constraints]
    >>>
    >>> # Solve the problem.
    >>> options = {'maxfev': 200}
    >>> res = minimize(lambda x: x[0]**2 + x[1]**2, [0, 1], bounds=bounds, constraints=constraints, options=options)
    >>> res.x
    array([0. , 0.5])
    '''

    linear_constraint, nonlinear_constraint_function = process_constraints(constraints)

    options = {'quiet': True} if options is None else options
    quiet = options.get("quiet", True)

    if method is None:
        if nonlinear_constraint_function is not None:
            if not quiet: print("Nonlinear constraints detected, applying COBYLA")
            method = "cobyla"
        elif linear_constraint is not None:
            if not quiet: print("Linear constraints detected without nonlinear constraints, applying LINCOA")
            method = "lincoa"
        elif bounds is not None:
            if not quiet: print("Bounds without linear or nonlinear constraints detected, applying BOBYQA")
            method = "bobyqa"
        else:
            if not quiet: print("No bounds or constraints detected, applying NEWUOA")
            method = "newuoa"
    else:
        # Raise some errors if methods were called with inappropriate options
        method = method.lower()
        if method not in ('newuoa', 'uobyqa', 'bobyqa', 'cobyla', 'lincoa'):
            raise ValueError(f"Method must be one of NEWUOA, UOBYQA, BOBYQA, COBYLA, or LINCOA, not '{method}'")
        if method != "cobyla" and nonlinear_constraint_function is not None:
            raise ValueError("Nonlinear constraints were provided for an algorithm that cannot handle them")
        if method not in ("cobyla", "lincoa") and linear_constraint is not None:
            raise ValueError("Linear constraints were provided for an algorithm that cannot handle them")
        if method not in ("cobyla", "bobyqa", "lincoa") and bounds is not None:
            raise ValueError("Bounds were provided for an algorithm that cannot handle them")

    try:
        lenx0 = len(x0)
        x0_is_scalar = False
    except TypeError:
        lenx0 = 1
        x0_is_scalar = True

    lb, ub = process_bounds(bounds, lenx0)

    # Check which variables are fixed and eliminate them from the problem.
    # Save the indices and values so that we can call the original function with
    # an array of the appropriate size, and so that we can add the fixed values to the
    # result when COBYLA returns.
    tol = get_arrays_tol(lb, ub)
    _fixed_idx = (
        (lb <= ub)
        & (ub <= lb + tol)
    )
    if any(_fixed_idx) and not all(_fixed_idx):
        # We should NOT reduce the problem if all variables are fixed. Otherwise, Aineq would be [], and
        # then bineq will be set to [] in the end. In this way, we lose completely the information in
        # these constraints. Consequently, we cannot evaluate the constraint violation correctly when needed.
        _fixed_values = 0.5 * (
            lb[_fixed_idx] + ub[_fixed_idx]
        )
        _fixed_values = np.clip(
            _fixed_values,
            lb[_fixed_idx],
            ub[_fixed_idx],
        )
        if isinstance(x0, np.ndarray):
            x0 = np.array(x0)[~_fixed_idx]
        else:
            # In this case x is presumably a list, so we turn it into a numpy array
            # for the convenience of indexing and then turn it back into a list
            x0 = np.array(x0)[~_fixed_idx].tolist()
        lb = lb[~_fixed_idx]
        ub = ub[~_fixed_idx]
        original_fun = fun
        def fixed_fun(x):
            newx = np.zeros(lenx0)
            newx[_fixed_idx] = _fixed_values
            newx[~_fixed_idx] = x
            return original_fun(newx, *args)
        fun = fixed_fun
        # Adjust linear_constraint
        if linear_constraint:
            new_lb = linear_constraint.lb - linear_constraint.A[:, _fixed_idx] @ _fixed_values
            new_ub = linear_constraint.ub - linear_constraint.A[:, _fixed_idx] @ _fixed_values
            new_A = linear_constraint.A[:, ~_fixed_idx]
            linear_constraint = LinearConstraint(new_A, new_lb, new_ub)
        # Adjust nonlinear constraints
        if nonlinear_constraint_function:
            original_nlc_fun = nonlinear_constraint_function
            def fixed_nlc_fun(x):
                newx = np.zeros(lenx0)
                newx[_fixed_idx] = _fixed_values
                newx[~_fixed_idx] = x
                return original_nlc_fun(newx, *args)
            nonlinear_constraint_function = fixed_nlc_fun


    # Project x0 onto the feasible set
    if nonlinear_constraint_function is None:
        result = _project(x0, lb, ub, {"linear": linear_constraint, "nonlinear": None})
        x0 = result.x
        # _project will upgrade x0 to a 1D array if it was a scalar, but the objective function
        # might expect a scalar, so we downgrade it back to a scalar if that's what it was originally
        if x0_is_scalar:
            x0 = x0[0]
    
    if linear_constraint is not None:
        A_eq, b_eq, A_ineq, b_ineq = separate_LC_into_eq_and_ineq(linear_constraint)
    else:
        A_eq, b_eq, A_ineq, b_ineq = None, None, None, None

    if all(_fixed_idx):
        x = 0.5 * (
            lb[_fixed_idx] + ub[_fixed_idx]
        )
        x = np.clip(
            x,
            lb[_fixed_idx],
            ub[_fixed_idx],
        )
        success = True
        status = FIXED_SUCCESS
        message = "All variables were fixed by the provided bounds."
        fun = fun(x)
        nfev = 1
        nlconstr = nonlinear_constraint_function(x) if nonlinear_constraint_function is not None else np.array([])
        maxcv = max(max((A_ineq @ x) - b_ineq) if A_ineq is not None else 0,
                    max((abs((A_eq @ x) - b_eq))) if A_eq is not None else 0,
                    max(np.append(0, nlconstr)))
        result = OptimizeResult(
            x = x,
            success = success,
            status = status,
            message = message,
            fun = fun,
            nfev = nfev,
            maxcv = maxcv,
            nlconstr = nlconstr,
            method = method
        )
        return result

    if nonlinear_constraint_function is not None:
        # If there is a nonlinear constraint function, we will call COBYLA, which needs the number of nonlinear
        # constraints (m_nlcon). In order to get this number we need to evaluate the constraint function at x0.
        # The constraint value at x0 (nlconstr0) is not discarded but passed down to the backend, as its
        # evaluation is assumed to be expensive. We also evaluate the objective function at x0 and pass the result
        # (f0) down to the backend, which expects nlconstr0 and f0 to be provided in sync.

        f0 = fun(x0, *args)
        nlconstr0 = nonlinear_constraint_function(x0)
        m_nlcon = len(nlconstr0)
    else:
        f0 = None
        nlconstr0 = None
        m_nlcon = 0

    if 'backend' not in options:  # default to Fortran
        options['backend'] = 'Fortran'

    if method.lower().strip() != 'cobyla' and options['backend'].lower() == 'python':
        warn('The pure Python implementation only supports COBYLA at this time. '
             'The Fortran implementation will be used instead.')
        options['backend'] = 'Fortran'

    if options['backend'].lower() == 'fortran':
        from .backends.pybindings import minimize as _minimize
        if m_nlcon > 0:
            options['f0'] = f0
            options['nlconstr0'] = nlconstr0
            options['m_nlcon'] = m_nlcon
        result = _minimize(
            fun,
            x0,
            args,
            method,
            lb,
            ub,
            A_eq,
            b_eq,
            A_ineq,
            b_ineq,
            nonlinear_constraint_function,
            callback,
            options
        )
        result = OptimizeResult(
            x = result.x,
            success = result.success,
            status = result.status,
            message = result.message,
            fun = result.fun,
            nfev = result.nfev,
            maxcv = result.maxcv,
            nlconstr = result.nlconstr,
            method = result.method,
        )
    elif options['backend'].lower() == "python":
        from .backends.pyprima.cobyla.cobyla import cobyla
        def calcfc(x):
            f = fun(x, *args)
            if nonlinear_constraint_function is not None:
                nlconstr = nonlinear_constraint_function(x)
            else:
                nlconstr = np.zeros(0)
            return f, nlconstr
        options.pop('backend', None)
        result = cobyla(
            calcfc,
            m_nlcon,
            x0,
            A_ineq,
            b_ineq,
            A_eq,
            b_eq,
            lb,
            ub,
            f0=f0,
            nlconstr0=nlconstr0,
            callback=callback,
            **options
        )
        result = OptimizeResult(
            x = result.x,
            success = result.info == 0,  # TODO: No magic numbers
            status = result.info,
            # message = result.message,
            fun = result.f,
            nfev = result.nf,
            maxcv = result.cstrv,
            nlconstr = result.constr,
            method = method,
        )
    else:
        raise ValueError(f"Backend must be either 'Fortran' or 'Python', not '{options['backend']}'")

    if any(_fixed_idx):
        newx = np.zeros(lenx0) + np.nan
        newx[_fixed_idx] = _fixed_values
        newx[~_fixed_idx] = result.x
        result.x = newx
    return result
