#ifndef PRIMA_CPP_COBYLA_GEOMETRY_HPP
#define PRIMA_CPP_COBYLA_GEOMETRY_HPP

// This module contains subroutines concerning the geometry-improving of the interpolation set.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <tuple>

#include "../consts.hpp"
#include "../linalg.hpp"

namespace prima {

inline int setdrop_tr(bool ximproved, const Eigen::VectorXd& d, double delta,
                      double rho, const Eigen::MatrixXd& sim, const Eigen::MatrixXd& simi) {
    // This function finds (the index) of a current interpolation point to be replaced with
    // the trust-region trial point. See (19)--(22) of the COBYLA paper.
    // N.B.:
    // 1. If XIMPROVED == True, then JDROP > 0 so that D is included into XPT. Otherwise,
    //    it is a bug.
    // 2. COBYLA never sets JDROP = NUM_VARS
    // TODO: Check whether it improves the performance if JDROP = NUM_VARS is allowed when
    // XIMPROVED is True. Note that UPDATEXFC should be revised accordingly.

    int num_vars = sim.rows();
    int jdrop = -1;

    // DISTSQ[j] is the square of the distance from the jth vertex of the simplex to get "best" point so
    // far, taking the trial point SIM[:, NUM_VARS] + D into account.
    Eigen::VectorXd distsq(sim.cols());
    if (ximproved) {
        for (int i = 0; i < num_vars; ++i)
            distsq[i] = (sim.col(i) - d).squaredNorm();
        distsq[num_vars] = d.squaredNorm();
    } else {
        for (int i = 0; i < num_vars; ++i)
            distsq[i] = sim.col(i).squaredNorm();
        distsq[num_vars] = 0;
    }

    // -------------------------------------------------------------------------------------------------- //
    //  The following code is Powell's scheme for defining JDROP.
    // -------------------------------------------------------------------------------------------------- //
    // ! JDROP = 0 by default. It cannot be removed, as JDROP may not be set below in some cases (e.g.,
    // ! when XIMPROVED == FALSE, MAXVAL(ABS(SIMID)) <= 1, and MAXVAL(VETA) <= EDGMAX).
    // jdrop = 0
    //
    // ! SIMID(J) is the value of the J-th Lagrange function at D. It is the counterpart of VLAG in UOBYQA
    // ! and DEN in NEWUOA/BOBYQA/LINCOA, but it excludes the value of the (N+1)-th Lagrange function.
    // simid = matprod(simi, d)
    // if (any(abs(simid) > 1) .or. (ximproved .and. any(.not. is_nan(simid)))) then
    //     jdrop = int(maxloc(abs(simid), mask=(.not. is_nan(simid)), dim=1), kind(jdrop))
    //     !!MATLAB: [~, jdrop] = max(simid, [], 'omitnan');
    // end if
    //
    // ! VETA(J) is the distance from the J-th vertex of the simplex to the best vertex, taking the trial
    // ! point SIM(:, N+1) + D into account.
    // if (ximproved) then
    //     veta = sqrt(sum((sim(:, 1:n) - spread(d, dim=2, ncopies=n))**2, dim=1))
    //     !!MATLAB: veta = sqrt(sum((sim(:, 1:n) - d).^2));  % d should be a column! Implicit expansion
    // else
    //     veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
    // end if
    //
    // ! VSIG(J) (J=1, .., N) is the Euclidean distance from vertex J to the opposite face of the simplex.
    // vsig = ONE / sqrt(sum(simi**2, dim=2))
    // sigbar = abs(simid) * vsig
    //
    // ! The following JDROP will overwrite the previous one if its premise holds.
    // mask = (veta > factor_delta * delta .and. (sigbar >= factor_alpha * delta .or. sigbar >= vsig))
    // if (any(mask)) then
    //     jdrop = int(maxloc(veta, mask=mask, dim=1), kind(jdrop))
    //     !!MATLAB: etamax = max(veta(mask)); jdrop = find(mask & ~(veta < etamax), 1, 'first');
    // end if
    //
    // ! Powell's code does not include the following instructions. With Powell's code, if SIMID consists
    // ! of only NaN, then JDROP can be 0 even when XIMPROVED == TRUE (i.e., D reduces the merit function).
    // ! With the following code, JDROP cannot be 0 when XIMPROVED == TRUE, unless VETA is all NaN, which
    // ! should not happen if X0 does not contain NaN, the trust-region/geometry steps never contain NaN,
    // ! and we exit once encountering an iterate containing Inf (due to overflow).
    // if (ximproved .and. jdrop <= 0) then  ! Write JDROP <= 0 instead of JDROP == 0 for robustness.
    //     jdrop = int(maxloc(veta, mask=(.not. is_nan(veta)), dim=1), kind(jdrop))
    //     !!MATLAB: [~, jdrop] = max(veta, [], 'omitnan');
    // end if
    // -------------------------------------------------------------------------------------------------- //
    //  Powell's scheme ends here.
    // -------------------------------------------------------------------------------------------------- //

    // The following definition of JDROP is inspired by SETDROP_TR in UOBYQA/NEWUOA/BOBYQA/LINCOA.
    // It is simpler and works better than Powell's scheme. Note that we allow JDROP to be NUM_VARS+1 if
    // XIMPROVED is True, whereas Powell's code does not.
    // See also (4.1) of Scheinberg-Toint-2010: Self-Correcting Geometry in Model-Based Algorithms for
    // Derivative-Free Unconstrained Optimization, which refers to the strategy here as the "combined
    // distance/poisedness criteria".

    double denom = std::max(rho, delta / 10.0);
    double denom_sq = denom * denom;
    Eigen::VectorXd weight(distsq.size());
    for (int i = 0; i < weight.size(); ++i)
        weight[i] = std::max(1.0, distsq[i] / denom_sq);
    // Similar to Powell's NEWUOA code.

    // Other possible definitions of weight. They work almost the same as the one above.
    // weight = distsq  // Similar to Powell's LINCOA code, but WRONG. See comments in LINCOA/geometry.f90.
    // weight = max(1, max(25 * distsq / delta**2))  // Similar to Powell's BOBYQA code, works well.
    // weight = max(1, max(10 * distsq / delta**2))
    // weight = max(1, max(1e2 * distsq / delta**2))
    // weight = max(1, max(distsq / rho**2))  ! Similar to Powell's UOBYQA

    // If 0 <= j < NUM_VARS, SIMID[j] is the value of the jth Lagrange function at D; the value of the
    // (NUM_VARS+1)th Lagrange function is 1 - sum(SIMID). [SIMID, 1 - sum(SIMID)] is the counterpart of
    // VLAG in UOBYQA and DEN in NEWUOA/BOBYQA/LINCOA.
    Eigen::VectorXd simid = matprod(simi, d);
    double simid_sum = simid.sum();
    Eigen::VectorXd score(sim.cols());
    for (int i = 0; i < num_vars; ++i)
        score[i] = weight[i] * std::abs(simid[i]);
    score[num_vars] = weight[num_vars] * std::abs(1.0 - simid_sum);

    // If XIMPROVED = False (D does not render a better X), set SCORE[NUM_VARS] = -1 to avoid JDROP = NUM_VARS.
    if (!ximproved)
        score[num_vars] = -1;

    // score[j] is NaN implies SIMID[j] is NaN, but we want abs(SIMID) to be big. So we
    // exclude such j.
    for (int i = 0; i < score.size(); ++i) {
        if (std::isnan(score[i]))
            score[i] = -1;
    }

    // The following if statement works a bit better than
    // `if any(score > 1) or (any(score > 0) and ximproved)` from Powell's UOBYQA and
    // NEWUOA code.
    bool any_positive = false;
    for (int i = 0; i < score.size(); ++i) {
        if (score[i] > 0) { any_positive = true; break; }
    }

    if (any_positive) {  // Powell's BOBYQA and LINCOA code.
        double max_score = -1;
        for (int i = 0; i < score.size(); ++i) {
            if (score[i] > max_score) {
                max_score = score[i];
                jdrop = i;
            }
        }
    }

    // JDROP >= 1 when XIMPROVED = TRUE unless NaN occurs in DISTSQ, which should not happen if the
    // starting point does not contain NaN and the trust-region/geometry steps never contain NaN.
    if (ximproved && jdrop < 0) {
        double max_distsq = -1;
        for (int i = 0; i < distsq.size(); ++i) {
            if (distsq[i] > max_distsq) {
                max_distsq = distsq[i];
                jdrop = i;
            }
        }
    }

    return jdrop;
}

inline Eigen::VectorXd geostep(int jdrop, const Eigen::MatrixXd* amat,
                                const Eigen::VectorXd* bvec,
                                const Eigen::MatrixXd& conmat, double cpen,
                                const Eigen::VectorXd& cval, double delbar,
                                const Eigen::VectorXd& fval, const Eigen::MatrixXd& simi) {
    // This function calculates a geometry step so that the geometry of the interpolation set is improved
    // when SIM[: JDROP_GEO] is replaced with SIM[:, NUM_VARS] + D. See (15)--(17) of the COBYLA paper.

    int m_lcon = (bvec != nullptr) ? bvec->size() : 0;
    int num_constraints = conmat.rows();
    int num_vars = simi.rows();

    // SIMI[JDROP, :] is a vector perpendicular to the face of the simplex to the opposite of vertex
    // JDROP. Set D to the vector in this direction and with length DELBAR.
    Eigen::VectorXd d = simi.row(jdrop).transpose();
    d = delbar * (d / norm(d));

    // The code below chooses the direction of D according to an approximation of the merit function.
    // See (17) of the COBYLA paper and line 225 of Powell's cobylb.f.

    // Calculate the coefficients of the linear approximations to the objective and constraint functions.
    // N.B.: CONMAT and SIMI have been updated after the last trust-region step, but G and A have not.
    // So we cannot pass G and A from outside.
    Eigen::VectorXd fdiff = (fval.head(num_vars).array() - fval[num_vars]).matrix();
    Eigen::VectorXd g = matprod(fdiff, simi);

    Eigen::MatrixXd A(num_vars, num_constraints);
    A.setZero();
    if (amat != nullptr && m_lcon > 0) {
        A.leftCols(m_lcon) = amat->transpose();
    }
    int m_nlcon = num_constraints - m_lcon;
    if (m_nlcon > 0) {
        Eigen::MatrixXd diff = conmat.block(m_lcon, 0, m_nlcon, num_vars).colwise()
                               - conmat.col(num_vars).segment(m_lcon, m_nlcon);
        A.rightCols(m_nlcon) = matprod(diff, simi).transpose();
    }

    // CVPD and CVND are the predicted constraint violation of D and -D by the linear models.
    double cvpd = 0;
    double cvnd = 0;
    for (int i = 0; i < num_constraints; ++i) {
        double adi = inprod(d, A.col(i));
        cvpd = std::max(cvpd, conmat(i, num_vars) + adi);
        cvnd = std::max(cvnd, conmat(i, num_vars) - adi);
    }

    if (-inprod(d, g) + cpen * cvnd < inprod(d, g) + cpen * cvpd)
        d *= -1;

    return d;
}

} // namespace prima

#endif // PRIMA_CPP_COBYLA_GEOMETRY_HPP
