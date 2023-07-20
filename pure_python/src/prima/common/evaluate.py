import numpy as np
from prima.common.consts import FUNCMAX, CONSTRMAX, REALMAX, DEBUGGING

# This is a module evaluating the objective/constraint function with
# Nan/Inf handling.


def moderatex(x):
    '''
    This function moderates a decision variable. It replaces NaN by 0 and Inf/-Inf by REALMAX/-REALMAX.
    '''
    y = x
    y[np.isnan(x)] = 0
    y = np.maximum(-REALMAX, np.minimum(REALMAX, y))
    return y

def moderatef(f):
    f = FUNCMAX if np.isnan(f) else f
    return min(FUNCMAX, f)

def moderatec(c):
    np.nan_to_num(c, copy=False, nan=-CONSTRMAX)
    c = np.clip(c, -CONSTRMAX, CONSTRMAX)
    return c


def evaluate(calcfc, x):
    """
    This function evaluates CALCFC at X, setting F to the objective function value, CONSTR to the
    constraint value, and CSTRV to the constraint violation. Nan/Inf are handled by a moderated
    extreme barrier.
    """

    # Preconditions
    if DEBUGGING:
        # X should not contain NaN if the initial X does not contain NaN and the subroutines generating
        # trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
        assert not any(np.isnan(x))

    #====================#
    # Calculation starts #
    #====================#

    if any(np.isnan(x)):
        f = np.sum(x)
        constr = f  # TODO: This is supposed to be an array, but I don't know how long it is here.
        cstrv = f
    else:
        f, constr = calcfc(moderatex(x))  # Evaluate F and CONSTR; We moderate X before doing so.
        
        # Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
        # finite value. This is naive, and better approaches surely exist.
        f = moderatef(f)
        constr = moderatec(constr)

        # Evaluate the constraint violation for constraints CONSTR(X) >= 0.
        cstrv = np.max(np.append(-constr, 0))

    #==================#
    # Calculation ends #
    #==================#

    # Postconditions
    if DEBUGGING:
        # With X not containing NaN, and with the moderated extreme barrier, F cannot be NaN/+Inf, and
        # CONSTR cannot be NaN/-Inf.
        assert not (np.isnan(f) or np.isposinf(f))
        assert not any(np.isnan(constr) | np.isneginf(constr))
        assert not (cstrv < 0 or np.isnan(cstrv) or np.isposinf(cstrv))

    return f, constr, cstrv