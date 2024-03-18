import numpy as np
        
def process_single_nl_constraint(nlc, x0):
    '''
    The Python interfaces receives the constraints as lb <= constraint(x) <= ub, 
    but the Fortran backend expects the nonlinear constraints to be constraint(x) <= 0.
    Thus a conversion is needed.
    
    In addition to the conversion, we would like to check the size of the output of the
    constraint function as compared to the provided lb and ub, and so we must run the
    constraint function here. Since we presume it is expensive to run the constraint function,
    we must also return the output of the constraint function to the caller, so that they may
    avoid running the constraint function unnecessarily.
    
    In order to accomplish all these goals, we first run the constraint function, then we
    upgrade lb/ub to vectors if they are not already (while checking their sizes in the process
    and raising exceptions if necessary), and then we create a function to transform the output
    of the constraint function to the form expected by the Fortran backend.

    '''
    
    # Run the constraint function and upgrade the lb/ub to the correct size,
    # while also checking that if they were already vectors that they had the correct
    # size in the first place.
    nlconstr0 = nlc.fun(x0)
    nlconstr0 = np.atleast_1d(np.array(nlconstr0, dtype=np.float64))
    lb, ub = _upgrade_lb_ub_to_vectors(nlc.lb, nlc.ub, nlconstr0)
    
    # Check if any bounds contain -inf and +inf simultaneously
    if np.any((lb == -np.inf) & (ub == np.inf)):
        raise ValueError("A NonlinearConstraint was provided without specifying lower or upper bounds")
    
    align_constraint_values = lambda values: np.concatenate(([vi - ub_ii for ub_ii, vi in zip(ub, values) if ub_ii < np.inf],
                                [lb_ii - vi for lb_ii, vi in zip(lb, values) if lb_ii > -np.inf]))
    def newconstraint(x):
        values = np.atleast_1d(np.array(nlc.fun(x), dtype=np.float64))
        return align_constraint_values(values)
    nlconstr0 = align_constraint_values(nlconstr0)
    return newconstraint, nlconstr0


def process_nl_constraints(nlcs, x0):
    functions = []
    nlconstr0 = np.empty(0)
    for nlc in nlcs:
        fun_i, nlconstr0_i = process_single_nl_constraint(nlc, x0)
        functions.append(fun_i)
        nlconstr0 = np.concatenate((nlconstr0, nlconstr0_i))
    def constraint_function(x):
        values = np.empty(0)
        for fun in functions:
            values = np.concatenate((values, fun(x)))
        return values
    return constraint_function, nlconstr0


def _upgrade_lb_ub_to_vectors(lb, ub, nlconstr0):
    '''
    Make sure length of lb/ub align with length of nlconstr0
    '''
    lb_is_scalar = not hasattr(lb, "__len__")
    ub_is_scalar = not hasattr(ub, "__len__")
    len_nlconstr0 = len(nlconstr0)
    if lb_is_scalar and ub_is_scalar:
        return np.array([lb]*len_nlconstr0, dtype=np.float64), np.array([ub]*len_nlconstr0, dtype=np.float64)
    elif lb_is_scalar and not ub_is_scalar:
        if len(ub) != len_nlconstr0:
            raise ValueError("The number of elements in the constraint function's output does not match the number of elements in the upper bound.")
        return np.array([lb]*len_nlconstr0, dtype=np.float64), np.array(ub, dtype=np.float64)
    elif not lb_is_scalar and ub_is_scalar:
        if len(lb) != len_nlconstr0:
            raise ValueError("The number of elements in the constraint function's output does not match the number of elements in the lower bound.")
        return np.array(lb, dtype=np.float64), np.array([ub]*len_nlconstr0, dtype=np.float64)
    else:
        if len(lb) != len_nlconstr0:
            raise ValueError("The number of elements in the constraint function's output does not match the number of elements in the lower bound.")
        if len(ub) != len_nlconstr0:
            raise ValueError("The number of elements in the constraint function's output does not match the number of elements in the upper bound.")
        return np.array(lb, dtype=np.float64), np.array(ub, dtype=np.float64)
