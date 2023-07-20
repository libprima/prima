from prima.common.evaluate import evaluate, moderatex
from prima.common.consts import EPS, RHOBEG_DEFAULT, RHOEND_DEFAULT, \
    CTOL_DEFAULT, CWEIGHT_DEFAULT, FTARGET_DEFAULT, IPRINT_DEFAULT, \
    MAXFUN_DIM_DEFAULT, DEBUGGING
from prima.common.preproc import preproc
from prima.common.present import present
from prima.cobyla.cobylb import cobylb
import numpy as np
from dataclasses import dataclass


# Return type for COBYLA
@dataclass
class COBYLAResult:
    x: np.ndarray
    f: float
    constr: np.ndarray
    cstrv: float
    nf: int
    xhist: np.ndarray | None
    fhist: np.ndarray | None
    chist: np.ndarray | None
    conhist: np.ndarray | None
    info: int


def cobyla(calcfc, num_constraints, x, f0=None, constr0=None, rhobeg=None,
           rhoend=None, ftarget=FTARGET_DEFAULT, ctol=CTOL_DEFAULT, cweight=CWEIGHT_DEFAULT,
           maxfun=None, iprint=IPRINT_DEFAULT, eta1=None, eta2=None, gamma1=0.5, gamma2=2,
           maxhist=None, maxfilt=2000):
    """
    Among all the arguments, only CALCFC, M, and X are obligatory. The others are OPTIONAL and you
    can neglect them unless you are familiar with the algorithm. Any unspecified optional input will
    take the default value detailed below. For instance, we may invoke the solver as follows.

    First define CALCFC, NUM_CONSTRAINTS, and X, and then do the following.
    result = cobyla(calcfc, num_constraints, x)

    or

    First define CALCFC, NUM_CONSTRAINTS, and X, and then do the following.
    result = cobyla(calcfc, num_constraints, x, rhobeg=0.5, rhoend=1E-3, maxfun=100)

    IMPORTANT NOTICE: The user must set NUM_CONSTRAINTS correctly to the number of constraints, namely
    the size of CONSTR introduced below. Set NUM_CONSTRAINTS to 0 if there is no constraint.

    See examples/cobyla_example.py for a concrete example.

    A detailed introduction to the arguments is as follows. The section below is largely similar to the
    comments in the Fortran version but has been slightly modified to reflect the differences in
    some of the arguments.

    CALCFC
    F, CONSTR = CALCFC(X) should evaluate the objective function and constraints at the given
    vector X; it should return the objective function value and the constraint value. It must be
    provided by the user, and its definition must conform to the following interface:
    !-------------------------------------------------------------------------!
        def calcfc(x):
            # ... calculations ...
            return f, constr
    !-------------------------------------------------------------------------!
    X, F, and CONSTR are numpy arrays. If there are no constraints, CONSTR is expected to be an
    empty numpy array.

    NUM_CONSTRAINTS
    NUM_CONSTRAINTS must be set to the number of constraints
    N.B.:
    1. Why don't we define NUM_CONSTRAINTS as optional and determine it from running CALCFC?
    Because we need to know it ahead of time to setup up various matrices and vectors. Since CALCFC
    might be expensive, we prefer to avoid running it unnecessarily.

    X
    X should be an N-dimensional vector that contains the starting point, N being the
    dimension of the problem.

    F0
    F0, if present, should be set to the objective function value of the starting X.

    CONSTR0
    CONSTR0, if present, should be set to the constraint value of the starting X; in addition,
    SIZE(CONSTR0) must be NUM_CONSTRAINTS, or the solver will abort.

    RHOBEG, RHOEND
    Default: RHOBEG = 1, RHOEND = 10^-6. RHOBEG and RHOEND must be set to
    the initial and final values of a trust-region radius, both being positive and RHOEND <= RHOBEG.
    Typically RHOBEG should be about one tenth of the greatest expected change to a variable, and
    RHOEND should indicate the accuracy that is required in the final values of the variables.

    FTARGET
    Default: -Inf.
    FTARGET is the target function value. The algorithm will terminate when a feasible point with a
    function value <= FTARGET is found.

    CTOL
    Default: machine epsilon.
    CTOL is the tolerance of constraint violation. Any X with MAXVAL(-CONSTR(X)) <= CTOL is
    considered feasible.
    N.B.: 1. CTOL is absolute, not relative. 2. CTOL is used only when selecting the returned X.
    It does not affect the iterations of the algorithm.

    CWEIGHT
    Default: CWEIGHT_DFT defined in common/consts.py.
    CWEIGHT is the weight that the constraint violation takes in the selection of the returned X.

    MAXFUN
    Default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in common/consts.py.
    MAXFUN is the maximal number of calls of CALCFC.

    IPRINT
    Default: 0.
    The value of IPRINT should be set to 0, 1, -1, 2, -2, 3, or -3, which controls how much
    information will be printed during the computation:
    0: there will be no printing;
    1: a message will be printed to the screen at the return, showing the best vector of variables
        found and its objective function value;
    2: in addition to 1, each new value of RHO is printed to the screen, with the best vector of
        variables so far and its objective function value; each new value of CPEN is also printed;
    3: in addition to 2, each function evaluation with its variables will be printed to the screen;
    -1, -2, -3: the same information as 1, 2, 3 will be printed, not to the screen but to a file
        named COBYLA_output.txt; the file will be created if it does not exist; the new output will
        be appended to the end of this file if it already exists.
    Note that IPRINT = +/-3 can be costly in terms of time and/or space.

    ETA1, ETA2, GAMMA1, GAMMA2
    Default: ETA1 = 0.1, ETA2 = 0.7, GAMMA1 = 0.5, and GAMMA2 = 2.
    ETA1, ETA2, GAMMA1, and GAMMA2 are parameters in the updating scheme of the trust-region radius
    detailed in the subroutine TRRAD in trustregion.py. Roughly speaking, the trust-region radius
    is contracted by a factor of GAMMA1 when the reduction ratio is below ETA1, and enlarged by a
    factor of GAMMA2 when the reduction ratio is above ETA2. It is required that 0 < ETA1 <= ETA2
    < 1 and 0 < GAMMA1 < 1 < GAMMA2. Normally, ETA1 <= 0.25. It is NOT advised to set ETA1 >= 0.5.

    XHIST, FHIST, CHIST, CONHIST, MAXHIST
    MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
    XHIST, if present, will output the history of iterates; FHIST, if present, will output the
    history function values; CHIST, if present, will output the history of constraint violations;
    CONHIST, if present, will output the history of constraint values; MAXHIST should be a
    nonnegative integer, and XHIST/FHIST/CHIST/CONHIST will output only the history of the last
    MAXHIST iterations. Therefore, MAXHIST= 0 means XHIST/FHIST/CONHIST/CHIST will output nothing,
    while setting MAXHIST = MAXFUN requests XHIST/FHIST/CHIST/CONHIST to output all the history.
    If XHIST is present, its size at exit will be (N, min(NF, MAXHIST)); if FHIST is present, its
    size at exit will be min(NF, MAXHIST); if CHIST is present, its size at exit will be
    min(NF, MAXHIST); if CONHIST is present, its size at exit will be (M, min(NF, MAXHIST)).

    IMPORTANT NOTICE:
    Setting MAXHIST to a large value can be costly in terms of memory for large problems.
    MAXHIST will be reset to a smaller value if the memory needed exceeds MAXHISTMEM defined in
    CONSTS_MOD (see consts.F90 under the directory named "common").
    Use *HIST with caution!!! (N.B.: the algorithm is NOT designed for large problems).

    MAXFILT
    Input, INTEGER(IK) scalar.
    MAXFILT is a nonnegative integer indicating the maximal length of the filter used for selecting
    the returned solution; default: MAXFILT_DFT (a value lower than MIN_MAXFILT is not recommended);
    see common/consts.F90 for the definitions of MAXFILT_DFT and MIN_MAXFILT.

    Output: A Python dataclass with the following fields
    X
    X will be set to the best vector of variables found.

    F
    F will be set to the objective function value at X.

    CSTRV
    Optional, only appears if CSTRV is True. CSTRV will be set to the constraint violation of X at
    exit, i.e., MAXVAL([-CONSTR(X), 0]).

    CONSTR
    Optional, only appears if CONSTR is True. CONSTR will be set to the constraint value of X at
    exit.

    NF
    NF will be set to the number of calls of CALCFC at exit.

    INFO
    INFO is the exit flag. It will be set to one of the following values defined in the module
    INFOS_MOD (see common/infos.f90):
    SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
    FTARGET_ACHIEVED: the target function value is reached;
    MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
    MAXTR_REACHED: the trust region iteration has been performed MAXTR times (MAXTR = 2*MAXFUN);
    NAN_INF_X: NaN or Inf occurs in X;
    DAMAGING_ROUNDING: rounding errors are becoming damaging.
    !--------------------------------------------------------------------------!
    The following case(s) should NEVER occur unless there is a bug.
    NAN_INF_F: the objective function returns NaN or +Inf;
    NAN_INF_MODEL: NaN or Inf occurs in the model;
    TRSUBP_FAILED: a trust region step failed to reduce the model
    !--------------------------------------------------------------------------!
    """

    # Local variables
    solver = 'COBYLA'

    # Preconditions
    if DEBUGGING:
        assert (f0 is None) == (constr0 is None), "f0 and constr0 must be both present or both absent"

    # Sizes
    num_vars = len(x)

    # Exit if constr0 is not consistent with num_constraints
    if present(constr0):
        assert np.size(constr0) == num_constraints


    # If the user provides the function & constraint value at X0, then use them.
    # Otherwise, set [F0, CONSTR0] = [F(X0), CONSTR(X0)], so that COBYLB only needs a single interface.
    if (f0 is None and constr0 is None) or not all(np.isfinite(x)):
        # Replace any NaN in X by ZERO and Inf/-Inf in X by REALMAX/-REALMAX.
        x = moderatex(x)
        f0, constr0, cstrv0 = evaluate(calcfc, x)
    
    if present(rhobeg):
        rhobeg = rhobeg
    elif present(rhoend) and np.isfinite(rhoend) and rhoend > 0:
        rhobeg = np.maximum(10*rhoend, RHOBEG_DEFAULT)
    else:
        rhobeg = RHOBEG_DEFAULT

    if present(rhoend):
        rhoend = rhoend
    elif rhobeg > 0:
        rhoend = np.maximum(EPS, np.minimum(0.1*rhobeg, RHOEND_DEFAULT))
    else:
        rhoend = RHOEND_DEFAULT

    maxfun = maxfun if present(maxfun) else MAXFUN_DIM_DEFAULT*num_vars

    if present(eta1):
        eta1 = eta1
    elif present(eta2) and 0 < eta2 < 1:
        eta1 = np.maximum(EPS, eta2/7)
    else:
        eta1 = 0.1

    if present(eta2):
        eta2 = eta2
    else:
        if 0 < eta1 < 1:
            eta2 = (eta1 + 2)/3
        else:
            eta2 = 0.7

    maxhist = maxhist if present(maxhist) else np.max([maxfun, num_vars+2, MAXFUN_DIM_DEFAULT*num_vars])

    # Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
    # Python: The underscores refer to `npt` and `x0`, neither of which are relevant here.
    iprint, maxfun, maxhist, ftarget, rhobeg, rhoend, _, maxfilt, ctol, cweight, \
        eta1, eta2, gamma1, gamma2, _ = preproc(solver, num_vars, iprint, maxfun, maxhist, \
                                                 ftarget, rhobeg, rhoend, num_constraints=num_constraints, ctol=ctol, cweight=cweight, \
                                                 eta1=eta1, eta2=eta2, gamma1=gamma1, gamma2=gamma2, maxfilt=maxfilt)

    
    # Further revise MAXHIST according to MAXHISTMEM, and allocate memory for the history.
    # In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
    # CHIST = NaN(1, MAXFUN), CONHIST = NaN(M, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
    # if they are requested; replace MAXFUN with 0 for the history that is not requested.
    # prehist(maxhist, num_vars, present(xhist), xhist_loc, present(fhist), fhist_loc, &
    #     & present(chist), chist_loc, m, present(conhist), conhist_loc)
    

    # call coblyb, which performs the real calculations
    x, f, constr, cstrv, nf, xhist, fhist, chist, conhist, info = cobylb(calcfc, iprint, maxfilt, maxfun, ctol, cweight, eta1, eta2,
      ftarget, gamma1, gamma2, rhobeg, rhoend, constr0, f0, x, 
      cstrv0, maxhist)
    
    return COBYLAResult(x, f, constr, cstrv, nf, xhist, fhist, chist, conhist, info)
