import numpy as np

from warnings import warn  # TODO: Need to determine the text of the warning and actually warn about it


class NonlinearConstraint(object):
    def __init__(self, fun, lb=-np.inf, ub=0):
        self.fun = fun
        self.lb = lb
        self.ub = ub


def process_nl_constraints(x0, nlcs, options):
    '''
    This function will accomplish several things:
    1. It will evaluate all linear constraints at x0 to build up nlconstr0
    2. In doing so, it will determine the number of constraints
    3. It will lump multiple constraints into a single constraint whose function
       will call of the constraints in turn, and whose lb/ub are guaranteed to be
       numpy arrays with a length equal to the number of constraints
    '''
    num_constraints = []
    nlconstr0 = []
    for nlc in nlcs:
        num_constraints_i, nlconstr0_i = get_num_constraints_from_constraint_function(x0, nlc)
        num_constraints.append(num_constraints_i)
        try:
            nlconstr0.extend(nlconstr0_i)
        except TypeError:
            nlconstr0.append(nlconstr0_i)
    assert len(nlconstr0) == sum(num_constraints), "There is a mismatch in the detected number of constraints and the length of the constraints returned"
    options['nlconstr0'] = nlconstr0

    
    # Upgrade any potentially scalar lb/ub to vectors
    lb = []
    ub = []
    for nlc, num_constraints_i in zip(nlcs, num_constraints):
        try:
            lb.extend(nlc.lb)
        except TypeError:
            lb.extend([nlc.lb]*num_constraints_i)

        try:
            ub.extend(nlc.ub)
        except TypeError:
            ub.extend([nlc.ub]*num_constraints_i)
    # Turn them into numpy arrays to enable arithmetic operations
    lb = np.array(lb)
    ub = np.array(ub)

    # Combine all the constraint functions into a single function
    def nlc_fun(x):
        constraints = []
        for nlc in nlcs:
            constraints_i = nlc.fun(x)
            try:
                constraints.extend(constraints_i)
            except TypeError:
                constraints.append(constraints_i)
        return np.array(constraints)
    
    return NonlinearConstraint(nlc_fun, lb=lb, ub=ub)

    
def get_num_constraints_from_constraint_function(x0, nlc):
    nlconstr0 = nlc.fun(x0)
    try:
        num_constraints = len(nlconstr0)
    except TypeError:
        num_constraints = 1
    return num_constraints, nlconstr0
