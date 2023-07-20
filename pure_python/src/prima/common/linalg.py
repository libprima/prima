import numpy as np
from prima.common.consts import DEBUGGING, EPS, REALMAX, REALMIN
from prima.common.present import present

def planerot(x):
    '''
    As in MATLAB, planerot(x) returns a 2x2 Givens matrix G for x in R2 so that Y=G@x has Y[1] = 0.
    Roughly speaking, G = np.array([[x[0]/R, x[1]/R], [-x[1]/R, x[0]/R]]), where R = np.linalg.norm(x).
    0. We need to take care of the possibilities of R=0, Inf, NaN, and over/underflow.
    1. The G defined above is continuous with respect to X except at 0. Following this definition,
    G = np.array([[np.sign(x[0]), 0], [0, np.sign(x[0])]]) if x[1] == 0,
    G = np.array([[0, np.sign(x[1])], [np.sign(x[1]), 0]]) if x[0] == 0
    Yet some implementations ignore the signs, leading to discontinuity and numberical instability.
    2. Difference from MATLAB: if x contains NaN of consists of only Inf, MATLAB returns a NaN matrix,
    but we return an identity matrix or a matrix of +/-np.sqrt(2). We intend to keep G always orthogonal.
    '''

    # Preconditions
    if DEBUGGING:
        assert len(x) == 2, "x must be a 2-vector"

    # ==================
    # Calculation starts
    # ==================

    # Define C = X(1) / R and S = X(2) / R with R = HYPOT(X(1), X(2)). Handle Inf/NaN, over/underflow.
    if (any(np.isnan(x))):
        # In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
        c = 1
        s = 0
    elif (all(np.isinf(x))):
        # In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
        c = 1 / np.sqrt(2) * np.sign(x[0])
        s = 1 / np.sqrt(2) * np.sign(x[1])
    elif (abs(x[0]) <= 0 and abs(x[1]) <= 0): # X(1) == 0 == X(2).
        c = 1
        s = 0
    elif (abs(x[1]) <= EPS * abs(x[0])):
        # N.B.:
        # 0. With <= instead of <, this case covers X(1) == 0 == X(2), which is treated above separately
        # to avoid the confusing SIGN(., 0) (see 1).
        # 1. SIGN(A, 0) = ABS(A) in Fortran but sign(0) = 0 in MATLAB, Python, Julia, and R#
        # 2. Taking SIGN(X(1)) into account ensures the continuity of G with respect to X except at 0.
        c = np.sign(x[0])
        s = 0
    elif (abs(x[0]) <= EPS * abs(x[1])):
        # N.B.: SIGN(A, X) = ABS(A) * sign of X /= A * sign of X # Therefore, it is WRONG to define G
        # as SIGN(RESHAPE([ZERO, -ONE, ONE, ZERO], [2, 2]), X(2)). This mistake was committed on
        # 20211206 and took a whole day to debug# NEVER use SIGN on arrays unless you are really sure.
        c = 0
        s = np.sign(x[1])
    else:
        # Here is the normal case. It implements the Givens rotation in a stable & continuous way as in:
        # Bindel, D., Demmel, J., Kahan, W., and Marques, O. (2002). On computing Givens rotations
        # reliably and efficiently. ACM Transactions on Mathematical Software (TOMS), 28(2), 206-238.
        # N.B.: 1. Modern compilers compute SQRT(REALMIN) and SQRT(REALMAX/2.1) at compilation time.
        # 2. The direct calculation without involving T and U seems to work better; use it if possible.
        if (all(np.logical_and(np.sqrt(REALMIN) < np.abs(x), np.abs(x) < np.sqrt(REALMAX / 2.1)))):
            # Do NOT use HYPOTENUSE here; the best implementation for one may be suboptimal for the other
            r = np.linalg.norm(x)
            c = x[0] / r
            s = x[1] / r
        elif (abs(x[0]) > abs(x[1])):
            t = x[1] / x[0]
            u = max(1, abs(t), np.sqrt(1 + t**2))  # MAXVAL: precaution against rounding error.
            u *= np.sign(x[0]) ##MATLAB: u = sign(x(1))*sqrt(1 + t**2)
            c = 1 / u
            s = t / u
        else:
            t = x[0] / x[1]
            u = max([1, abs(t), np.sqrt(1 + t**2)])  # MAXVAL: precaution against rounding error.
            u *= np.sign(x[1]) ##MATLAB: u = sign(x(2))*sqrt(1 + t**2)
            c = t / u
            s = 1 / u

    G = np.array([[c, s], [-s, c]]) #  MATLAB: G = [c, s; -s, c]

    #====================#
    #  Calculation ends  #
    #====================#

    # Postconditions
    if DEBUGGING:
        assert G.shape == (2,2)
        assert np.all(np.isfinite(G))
        assert abs(G[0, 0] - G[1, 1]) + abs(G[0, 1] + G[1, 0]) <= 0
        tol = np.maximum(1.0E-10, np.minimum(1.0E-1, 1.0E6 * EPS))
        assert isorth(G, tol)
        if all(np.logical_and(np.isfinite(x), np.abs(x) < np.sqrt(REALMAX / 2.1))):
            r = np.linalg.norm(x)
            assert max(abs(G@x - [r, 0])) <= max(tol, tol * r), 'G @ X = [||X||, 0]'

    return G


def isminor(x, ref):
    '''
    This function tests whether x is minor compared to ref. It is used by Powell, e.g., in COBYLA.
    In precise arithmetic, isminor(x, ref) is true if and only if x == 0; in floating point
    arithmetic, isminor(x, ref) is true if x is 0 or its nonzero value can be attributed to
    computer rounding errrors according to ref.
    Larger sensitivity means the function is more strict/precise, the value 0.1 being due to Powell.

    For example:
    isminor(1e-20, 1e300) -> True, because in floating point arithmetic 1e-20 cannot be added to
    1e300 without being rounded to 1e300.
    isminor(1e300, 1e-20) -> False, because in floating point arithmetic adding 1e300 to 1e-20
    dominates the latter number.
    isminor(3, 4) -> False, because 3 can be added to 4 without being rounded off
    '''

    sensitivity = 0.1
    refa = abs(ref) + sensitivity * abs(x)
    refb = abs(ref) + 2 * sensitivity * abs(x)
    return np.logical_or(abs(ref) >= refa, refa >= refb)


def isinv(A, B, tol=None):
    '''
    This procedure tests whether A = B^{-1} up to the tolerance TOL.
    '''

    # Sizes
    n = np.size(A, 0)

    # Preconditions
    if DEBUGGING:
        assert np.size(A, 0) == np.size(A, 1)
        assert np.size(B, 0) == np.size(B, 1)
        assert np.size(A, 0) == np.size(B, 0)
        if present(tol):
            assert tol >= 0

    #====================#
    # Calculation starts #
    #====================#

    tol = tol if present(tol) else np.minimum(1e-3, 1e2 * EPS * np.maximum(np.size(A, 0), np.size(A, 1)))
    tol = np.max([tol, tol * np.max(abs(A)), tol * np.max(abs(B))])
    is_inv = ((abs(A@B) - np.eye(n)) <= tol).all() or ((abs(B@A - np.eye(n))) <= tol).all()

    #===================#
    #  Calculation ends #
    #===================#
    return is_inv


def isorth(A, tol=None):
    '''
    This function tests whether the matrix A has orthonormal columns up to the tolerance TOL.
    '''

    # Preconditions
    if DEBUGGING:
        if present(tol):
            assert tol >= 0

    #====================#
    # Calculation starts #
    #====================#

    num_vars = np.size(A, 1)

    if num_vars > np.size(A, 0):
        is_orth = False
    elif (np.isnan(np.sum(abs(A)))):
        is_orth = False
    else:
        if present(tol):
            is_orth = (abs(A.T@A - np.eye(num_vars)) <= np.maximum(tol, tol * np.max(abs(A)))).all()
        else:
            is_orth = (abs(A.T@A - np.eye(num_vars)) <= 0).all()

    #====================#
    #  Calculation ends  #
    #====================#
    return is_orth
