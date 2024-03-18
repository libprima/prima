import numpy as np
        
def process_single_nl_constraint(nlc):
    '''
    The Python interfaces receives the constraints as lb <= constraint(x) <= ub, 
    but the Fortran backend expects the nonlinear constraints to be constr(x) <= 0.
    Thus a conversion is needed.
    
    This function will take a NonlinearConstraint and return a new function which incorporates
    the lower and upper bounds in such a way as to satisfy the expectation of the Fortran
    backend. This new function first calls the original function from the provided NonlinearConstraint
    and obtains the values, and then passes those values to a second function which combines
    them with the lower and upper bounds as appropriate. The main work here is creating that
    second function.
    
    In general to convert lb <= constraint(x) <= ub to newconstraint(x) <= 0, all we need to do
    is define newconstraint(x) as [constraint(x) - ub, lb - constraint(x)]. There are some details
    regarding the type and content of lb and ub as detailed below.
    
    
    The upper and lower bounds could be either scalars or vectors, and so we have 4 possible
    cases to take into account:
    
    1. Both lb and ub are scalars
        In this case we have a further 4 cases!
        
        1a. lb == -inf and ub == inf
            This case makes no sense since it means that there are effectively no constraints,
            so we raise an exception.
        1b. lb == -inf and ub != inf
            This is a one sided constraint, so we can define newconstraint(x) as constraint(x) - ub
        1c. lb != -inf and ub == inf
            This is a one sided constraint, so we can define newconstraint(x) as lb - constraint(x)
        1d. lb != -inf and ub != inf
            Since we have both constraints we define newconstraint(x) as [constraint(x) - ub, lb - constraint(x)]
    2. lb is a scalar and ub is a vector
        In this case we have a further 2 cases.
        
        2a. lb == -inf
            This is a one sided constraint, so we can define newconstraint(x) as [constraint(x) - ub], however we
            should first check if there is any inf in ub and raise an exception since this would be same case as 1a.
        2b. lb != -inf
            In this case we can have inf in ub, so newconstraint needs to check for that and omit it from the constraints.
            See the code for the exact method of accomplishing this.
    3. lb is a vector and ub is a scalar
        This is the same as case 2, but with lb and ub reversed.
    4. Both lb and ub are vectors
        There are no subcases here, other than checking for -np.inf in lb and np.inf in ub and removing those constraints.
        However we also check for any particular index where ub is inf and lb is -inf simultaneously and raise an exception
        if so, since again this is like case 1a.

    '''
    lb_is_scalar = not hasattr(nlc.lb, "__len__")
    ub_is_scalar = not hasattr(nlc.ub, "__len__")
    if lb_is_scalar and ub_is_scalar:
        if nlc.lb == -np.inf and nlc.ub == np.inf:
            raise ValueError("A NonlinearConstraint was provided without specifying lower or upper bounds")
        elif nlc.lb == -np.inf:
            align_constraint_values = lambda values: values - nlc.ub
        elif nlc.ub == np.inf:
            align_constraint_values = lambda values: nlc.lb - values
        else:
            align_constraint_values = lambda values: np.concatenate((values - nlc.ub, nlc.lb - values))
    elif lb_is_scalar and not ub_is_scalar:
        if nlc.lb == -np.inf:
            if np.inf in nlc.ub:
                raise ValueError("A NonlinearConstraint was provided without specifying lower or upper bounds")
            align_constraint_values = lambda values: values - nlc.ub
        else:
            align_constraint_values = lambda values: np.concatenate(([vi - ub_ii for ub_ii, vi in zip(nlc.ub, values) if ub_ii < np.inf],
                                        nlc.lb - values))
    elif not lb_is_scalar and ub_is_scalar:
        if nlc.ub == np.inf:
            if -np.inf in nlc.lb:
                raise ValueError("A NonlinearConstraint was provided without specifying lower or upper bounds")
            align_constraint_values = lambda values: nlc.lb - values
        else:
            align_constraint_values = lambda values: np.concatenate((values - nlc.ub,
                    [lb_ii - vi for lb_ii, vi in zip(nlc.lb, values) if lb_ii > -np.inf]))
    else:
        if np.any((nlc.lb == -np.inf) & (nlc.ub == np.inf)):
            raise ValueError("A NonlinearConstraint was provided without specifying lower or upper bounds")
        align_constraint_values = lambda values: np.concatenate(([vi - ub_ii for ub_ii, vi in zip(nlc.ub, values) if ub_ii < np.inf],
                                   [lb_ii - vi for lb_ii, vi in zip(nlc.lb, values) if lb_ii > -np.inf]))
    def newconstraint(x):
        values = np.atleast_1d(np.array(nlc.fun(x), dtype=np.float64))
        return align_constraint_values(values)
    return newconstraint


def process_nl_constraints(nlcs):
    functions = []
    for nlc in nlcs:
        functions.append(process_single_nl_constraint(nlc))
    def constraint_function(x):
        values = np.empty(0)
        for fun in functions:
            values = np.concatenate((values, fun(x)))
        return values
    return constraint_function
