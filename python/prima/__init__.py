from ._prima import minimize as _minimize, __version__, PRIMAMessage
# Bounds may appear unused in this file but we need to import it to make it available to the user
from scipy.optimize import NonlinearConstraint, LinearConstraint, Bounds
from ._nonlinear_constraints import process_nl_constraints
from ._linear_constraints import (
    combine_multiple_linear_constraints,
    separate_LC_into_eq_and_ineq,
)
from ._bounds import process_bounds
from enum import Enum
import numpy as np
from ._common import _project


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


def process_constraints(constraints, x0):
    # First throw it back if it's an empty tuple
    if not constraints:
        return None, None, None
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
        nonlinear_constraint_function, nlconstr0 = process_nl_constraints(nonlinear_constraints, x0)
    else:
        nonlinear_constraint_function = None
        nlconstr0 = None

    # Determine if we have multiple linear constraints, just 1, or none, and process accordingly
    if len(linear_constraints) > 1:
        linear_constraint = combine_multiple_linear_constraints(linear_constraints)
    elif len(linear_constraints) == 1:
        linear_constraint = linear_constraints[0]
    else:
        linear_constraint = None

    return linear_constraint, nonlinear_constraint_function, nlconstr0


def minimize(fun, x0, args=(), method=None, bounds=None, constraints=(), callback=None, options=None):

    linear_constraint, nonlinear_constraint_function, nlconstr0 = process_constraints(constraints, x0)
        
    quiet = options.get("quiet", True) if options is not None else True

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
        
    if nonlinear_constraint_function is not None:
        m_nlcon = len(nlconstr0)
        f0 = fun(x0)
    else:
        m_nlcon = 0
        f0 = None

    return _minimize(
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
        options,
        nlconstr0,
        m_nlcon,
        f0,
    )
