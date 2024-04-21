from prima.common.consts import DEBUGGING
from prima.common.linalg import isinv
import numpy as np


def assess_geo(delta, factor_alpha, factor_beta, sim, simi):
    '''
    This function checks if an interpolation set has acceptable geometry as (14) of the COBYLA paper
    '''

    # Local variables
    itol = 0.1

    # Sizes
    num_vars = np.size(sim, 0)

    # Preconditions
    if DEBUGGING:
        assert num_vars >= 1
        assert np.size(sim, 0) == num_vars and np.size(sim, 1) == num_vars + 1
        assert np.isfinite(sim).all()
        assert all(np.max(abs(sim[:, :num_vars]), axis=0) > 0)
        assert np.size(simi, 0) == num_vars and np.size(simi, 1) == num_vars
        assert np.isfinite(simi).all()
        assert isinv(sim[:, :num_vars], simi, itol)
        assert delta > 0
        assert factor_alpha > 0 and factor_alpha < 1
        assert factor_beta > 1

    #====================#
    # Calculation starts #
    #====================#

    # Calculate the values of sigma and eta
    # veta[j] (0 <= j < num_vars) is the distance between the vertices j and 0 (the best vertex)
    # of the simplex.
    # vsig[j] (0 <= j < num_vars) is the distance from vertex j to its opposite face of the simplex.
    # Thus vsig <= veta.
    # N.B.: What about the distance from vertex N+1 to its opposite face? Consider the simplex
    # {V_{N+1}, V_{N+1} + L*e_1,... v_{N+1} + L*e_N}, where V_{N+1} is vertex N+1,
    # namely the current "best" point, [e_1, ..., e_n] is an orthogonal matrix, and L is a
    # constant in the order of delta. This simplex is optimal in the sense that the interpolation
    # system has the minimal condition number, i.e. 1. For this simplex, the distance from
    # V_{N+1} to its opposite face is L/sqrt(n).
    vsig = 1/np.sqrt(np.sum(simi**2, axis=1))
    veta = np.sqrt(np.sum(sim[:, :num_vars]**2, axis=0))
    adequate_geo = all(vsig >= factor_alpha * delta) and all(veta <= factor_beta * delta)

    #==================#
    # Calculation ends #
    #==================#

    return adequate_geo


def setdrop_tr(ximproved, d, delta, rho, sim, simi):
    '''
    This function finds (the index) of a current interpolation point to be replaced with the
    trust-region trial point. See (19)-(22) of the COBYLA paper.
    N.B.:
    1. If XIMPROVED == True, then JDROP >= 0 so that D is included into XPT. Otherwise, it is a bug.
    2. COBYLA never sets JDROP = NUM_VARS
    TODO: Check whether it improves the performance if JDROP = NUM_VARS is allowed when XIMPROVED is True.
    Note that UPDATEXFC should be revised accordingly.
    '''

    # Local variables
    itol = 0.1

    # Sizes
    num_vars = np.size(sim, 0)

    # Preconditions
    if DEBUGGING:
        assert num_vars >= 1
        assert np.size(d) == num_vars and all(np.isfinite(d))
        assert delta >= rho and rho > 0
        assert np.size(sim, 0) == num_vars and np.size(sim, 1) == num_vars + 1
        assert np.isfinite(sim).all()
        assert all(np.max(abs(sim[:, :num_vars]), axis=0) > 0)
        assert np.size(simi, 0) == num_vars and np.size(simi, 1) == num_vars
        assert np.isfinite(simi).all()
        assert isinv(sim[:, :num_vars], simi, itol)

    #====================#
    # Calculation starts #
    #====================#

    # -------------------------------------------------------------------------------------------------- #
    #  The following code is Powell's scheme for defining JDROP.
    # -------------------------------------------------------------------------------------------------- #
    # ! JDROP = 0 by default. It cannot be removed, as JDROP may not be set below in some cases (e.g.,
    # ! when XIMPROVED == FALSE, MAXVAL(ABS(SIMID)) <= 1, and MAXVAL(VETA) <= EDGMAX).
    # jdrop = 0
    # 
    # ! SIMID(J) is the value of the J-th Lagrange function at D. It is the counterpart of VLAG in UOBYQA
    # ! and DEN in NEWUOA/BOBYQA/LINCOA, but it excludes the value of the (N+1)-th Lagrange function.
    # simid = matprod(simi, d)
    # if (any(abs(simid) > 1) .or. (ximproved .and. any(.not. is_nan(simid)))) then
    #     jdrop = int(maxloc(abs(simid), mask=(.not. is_nan(simid)), dim=1), kind(jdrop))
    #     !!MATLAB: [~, jdrop] = max(simid, [], 'omitnan');
    # end if
    # 
    # ! VETA(J) is the distance from the J-th vertex of the simplex to the best vertex, taking the trial
    # ! point SIM(:, N+1) + D into account.
    # if (ximproved) then
    #     veta = sqrt(sum((sim(:, 1:n) - spread(d, dim=2, ncopies=n))**2, dim=1))
    #     !!MATLAB: veta = sqrt(sum((sim(:, 1:n) - d).^2));  % d should be a column! Implicit expansion
    # else
    #     veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
    # end if
    # 
    # ! VSIG(J) (J=1, .., N) is the Euclidean distance from vertex J to the opposite face of the simplex.
    # vsig = ONE / sqrt(sum(simi**2, dim=2))
    # sigbar = abs(simid) * vsig
    # 
    # ! The following JDROP will overwrite the previous one if its premise holds.
    # mask = (veta > factor_delta * delta .and. (sigbar >= factor_alpha * delta .or. sigbar >= vsig))
    # if (any(mask)) then
    #     jdrop = int(maxloc(veta, mask=mask, dim=1), kind(jdrop))
    #     !!MATLAB: etamax = max(veta(mask)); jdrop = find(mask & ~(veta < etamax), 1, 'first');
    # end if
    # 
    # ! Powell's code does not include the following instructions. With Powell's code, if SIMID consists
    # ! of only NaN, then JDROP can be 0 even when XIMPROVED == TRUE (i.e., D reduces the merit function).
    # ! With the following code, JDROP cannot be 0 when XIMPROVED == TRUE, unless VETA is all NaN, which
    # ! should not happen if X0 does not contain NaN, the trust-region/geometry steps never contain NaN,
    # ! and we exit once encountering an iterate containing Inf (due to overflow).
    # if (ximproved .and. jdrop <= 0) then  ! Write JDROP <= 0 instead of JDROP == 0 for robustness.
    #     jdrop = int(maxloc(veta, mask=(.not. is_nan(veta)), dim=1), kind(jdrop))
    #     !!MATLAB: [~, jdrop] = max(veta, [], 'omitnan');
    # end if
    # -------------------------------------------------------------------------------------------------- #
    #  Powell's scheme ends here.
    # -------------------------------------------------------------------------------------------------- #

    # The following definition of JDROP is inspired by SETDROP_TR in UOBYQA/NEWUOA/BOBYQA/LINCOA.
    # It is simpler and works better than Powell's scheme. Note that we allow JDROP to be NUM_VARS+1 if
    # XIMPROVED is True, whereas Powell's code does not.
    # See also (4.1) of Scheinberg-Toint-2010: Self-Correcting Geometry in Model-Based Algorithms for
    # Derivative-Free Unconstrained Optimization, which refers to the strategy here as the "combined
    # distance/poisedness criteria".

    # DISTSQ[j] is the square of the distance from the jth vertex of the simplex to get "best" point so
    # far, taking the trial point SIM[:, NUM_VARS] + D into account.
    distsq = np.zeros(np.size(sim, 1))
    if ximproved:
        distsq[:num_vars] = np.sum((sim[:, :num_vars] - np.tile(d, (num_vars, 1)).T)**2, axis=0)
        distsq[num_vars] = np.sum(d**2)
    else:
        distsq[:num_vars] = np.sum(sim[:, :num_vars]**2, axis=0)
        distsq[num_vars] = 0

    weight = np.maximum(1, distsq / np.maximum(rho, delta/10)**2)  # Similar to Powell's NEWUOA code.

    # Other possible definitions of weight. They work almost the same as the one above.
    # weight = distsq  # Similar to Powell's LINCOA code, but WRONG. See comments in LINCOA/geometry.f90.
    # weight = max(1, max(25 * distsq / delta**2))  # Similar to Powell's BOBYQA code, works well.
    # weight = max(1, max(10 * distsq / delta**2))
    # weight = max(1, max(1e2 * distsq / delta**2))
    # weight = max(1, max(distsq / rho**2))  ! Similar to Powell's UOBYQA

    # If 0 <= j < NUM_VARS, SIMID[j] is the value of the jth Lagrange function at D; the value of the
    # (NUM_VARS+1)th Lagrange function is 1 - sum(SIMID). [SIMID, 1 - sum(SIMID)] is the counterpart of
    # VLAG in UOBYQA and DEN in NEWUOA/BOBYQA/LINCOA.
    simid = simi@d
    score = weight * abs(np.array([*simid, 1 - np.sum(simid)]))

    # If XIMPORVED = False (D does not render a better X), set SCORE[NUM_VARS] = -1 to avoid JDROP = NUM_VARS.
    if not ximproved:
        score[num_vars] = -1
    
    # The following if statement works a bit better than `if any(score > 1) or (any(score > 0) and ximproved)`
    # from Powell's UOBYQA and NEWUOA code.
    if any(score > 0):  # Powell's BOBYQA and LINCOA code.
        jdrop = np.where(score == max(score[~np.isnan(score)]))[0][0]
    elif ximproved:
        jdrop = np.argmax(distsq)
    else:
        jdrop = None  # We arrive here when XIMPROVED = False and no entry of score is positive.

    #==================#
    # Calculation ends #
    #==================#

    # Postconditions
    if DEBUGGING:
        assert jdrop is None or (jdrop >= 0 and jdrop < num_vars + 1)
        assert jdrop <= num_vars or ximproved
        assert jdrop >= 0 or not ximproved
        # JDROP >= 1 when XIMPROVED = TRUE unless NaN occurs in DISTSQ, which should not happen if the
        # starting point does not contain NaN and the trust-region/geometry steps never contain NaN.

    return jdrop


def setdrop_geo(delta, factor_alpha, factor_beta, sim, simi):
    '''
    This function finds (the index) of a current interpolation point to be replaced with a
    geometry-improving point. See (15)-(16) of the COBYLA paper.
    N.B.: COBYLA never sets jdrop to NUM_VARS
    '''

    # Local variables
    itol = 0.1

    # Sizes
    num_vars = np.size(sim, 0)

    # Preconditions
    if DEBUGGING:
        assert num_vars >= 1
        assert np.size(sim, 0) == num_vars and np.size(sim, 1) == num_vars + 1
        assert np.isfinite(sim).all()
        assert all(np.max(abs(sim[:, :num_vars]), axis=0) > 0)
        assert np.size(simi, 0) == num_vars and np.size(simi, 1) == num_vars
        assert np.isfinite(simi).all()
        assert isinv(sim[:, :num_vars], simi, itol)
        assert factor_alpha > 0 and factor_alpha < 1
        assert factor_beta > 1
        assert not assess_geo(delta, factor_alpha, factor_beta, sim, simi)

    #====================#
    # Calculation starts #
    #====================#

    # Calculate the values of sigma and eta
    # VSIG[j] for j = 0...NUM_VARS-1 is the Euclidean distance from vertex J to the opposite face of the simplex.
    vsig = 1 / np.sqrt(np.sum(simi**2, axis=1))
    veta = np.sqrt(np.sum(sim[:, :num_vars]**2, axis=0))

    # Decide which vertex to drop from the simplex. It will be replaced with a new point to improve the
    # acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
    if any(veta > factor_beta * delta):
        jdrop = np.where(veta == max(veta[~np.isnan(veta)]))[0][0]
    elif any(vsig < factor_alpha * delta):
        jdrop = np.where(vsig == min(vsig[~np.isnan(vsig)]))[0][0]
    else:
        # We arrive here if vsig and veta are all nan, which can happen due to nan in sim and simi
        # which should not happen unless there is a bug
        jdrop = None

    # Zaikun 230202: What if we consider veta and vsig together? The following attempts do not work well.
    # jdrop = max(np.sum(sim[:, :num_vars]**2, axis=0)*np.sum(simi**2, axis=1))  # Condition number
    # jdrop = max(np.sum(sim[:, :num_vars]**2, axis=0)**2 * np.sum(simi**2, axis=1))  # Condition number times distance

    #==================#
    # Calculation ends #
    #==================#

    # Postconditions
    if DEBUGGING:
        assert jdrop >= 0 and jdrop < num_vars
    return jdrop


def geostep(jdrop, cpen, conmat, cval, delta, fval, factor_gamma, simi):
    '''
    This function calculates a geometry step so that the geometry of the interpolation set is improved
    when SIM[: JDROP_GEO] is replaced with SIM[:NUM_VARS] + D. See (15)--(17) of the COBYLA paper.
    '''

    # Sizes
    num_constraints = np.size(conmat, 0)
    num_vars = np.size(simi, 0)

    # Preconditions
    if DEBUGGING:
        assert num_constraints >= 0
        assert num_vars >= 1
        assert delta > 0
        assert cpen > 0
        assert np.size(simi, 0) == num_vars and np.size(simi, 1) == num_vars
        assert np.isfinite(simi).all()
        assert np.size(fval) == num_vars + 1 and not any(np.isnan(fval) | np.isposinf(fval))
        assert np.size(conmat, 0) == num_constraints and np.size(conmat, 1) == num_vars + 1
        assert not (np.isnan(conmat) | np.isneginf(conmat)).any()
        assert np.size(cval) == num_vars + 1 and not any(cval < 0 | np.isnan(cval) | np.isposinf(cval))
        assert jdrop >= 0 and jdrop < num_vars
        assert factor_gamma > 0 and factor_gamma < 1

    #====================#
    # Calculation starts #
    #====================#

    # SIMI[JDROP, :] is a vector perpendicular to the face of the simplex to the opposite of vertex
    # JDROP. Thus VSIGJ * SIMI[JDROP, :] is the unit vector in this direction
    vsigj = 1 / np.sqrt(np.sum(simi[jdrop, :]**2))

    # Set D to the vector in the above-mentioned direction and with length FACTOR_GAMMA * DELTA. Since
    # FACTOR_ALPHA < FACTOR_GAMMA < FACTOR_BETA, D improves the geometry of the simplex as per (14) of
    # the COBYLA paper. This also explains why this subroutine does not replace DELTA with
    # DELBAR = max(min(0.1 * np.sqrt(max(DISTSQ)), 0.5 * DELTA), RHO) as in NEWUOA.
    d = factor_gamma * delta * (vsigj * simi[jdrop, :])

    # Calculate the coefficients of the linear approximations to the objective and constraint functions,
    # placing minus the objective function gradient after the constraint gradients in the array A
    A = np.zeros((num_vars, num_constraints + 1))
    A[:, :num_constraints] = ((conmat[:, :num_vars] - np.tile(conmat[:, num_vars], (num_vars, 1)).T)@simi).T
    A[:, num_constraints] = (fval[num_vars] - fval[:num_vars])@simi
    cvmaxp = np.max(np.append(0, -d@A[:, :num_constraints] - conmat[:, num_vars]))
    cvmaxn = np.max(np.append(0, d@A[:, :num_constraints] - conmat[:, num_vars]))
    if 2 * np.dot(d, A[:, num_constraints]) < cpen * (cvmaxp - cvmaxn):
        d *= -1

    #==================#
    # Calculation ends #
    #==================#

    # Postconditions
    if DEBUGGING:
        assert np.size(d) == num_vars and all(np.isfinite(d))
        # In theory, ||S|| == FACTOR_GAMMA*DELTA, which may be false due to rounding, but not too far.
        # It is crucial to ensure that the geometry step is nonzero, which holds in theory.
        assert np.linalg.norm(d) > 0.9 * factor_gamma * delta and np.linalg.norm(d) <= 1.1 * factor_gamma * delta
    return d