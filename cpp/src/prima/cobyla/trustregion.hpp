#ifndef PRIMA_CPP_COBYLA_TRUSTREGION_HPP
#define PRIMA_CPP_COBYLA_TRUSTREGION_HPP

// This module provides subroutines concerning the trust-region calculations of COBYLA.
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
#include "../powalg.hpp"

namespace prima {

namespace detail {

// TRSTLP_SUB: the real calculations for trstlp, both stage 1 and stage 2.
// Major differences between stage 1 and stage 2:
// 1. Initialization. Stage 2 inherits the values of some variables from stage 1, so
//    they are initialized in stage 1 but not in stage 2.
// 2. cviol. Updated after each iteration in stage 1, constant in stage 2.
// 3. sdirn. See the definition in the code for details.
// 4. optnew. The two stages have different objectives.
// 5. step <= cviol in stage 1.
//
// Powell's code can encounter infinite cycling, which did happen when testing CUTEst
// problems: DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
// Indeed, in all these cases, Inf/NaN appear in d due to extremely large values in
// A (up to 10^219). To resolve this, we set the maximal number of iterations to
// maxiter, and terminate if Inf/NaN occurs in d.
inline std::tuple<Eigen::VectorXi, int, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>
trstlp_sub(Eigen::VectorXi iact, int nact, int stage,
           const Eigen::MatrixXd& A, const Eigen::VectorXd& b, double delta,
           const Eigen::VectorXd& d_in, const Eigen::VectorXd& vmultc_in,
           const Eigen::MatrixXd& z_in) {
    int mcon = A.cols();
    int num_vars = A.rows();

    Eigen::VectorXd zdasav = Eigen::VectorXd::Zero(z_in.cols());
    Eigen::VectorXd vmultd = Eigen::VectorXd::Zero(vmultc_in.size());
    Eigen::VectorXd zdota = Eigen::VectorXd::Zero(z_in.cols());

    Eigen::MatrixXd z = z_in;
    Eigen::VectorXd vmultc = vmultc_in;
    Eigen::VectorXd d = d_in;

    int num_constraints = 0;
    int icon = 0;
    double cviol = 0;

    // Initialize according to stage
    if (stage == 1) {
        for (int i = 0; i < mcon; ++i) iact[i] = i;
        nact = 0;
        d.setZero();
        if (mcon == 0) {
            // Quick return: no constraints in stage 1
            return {iact, nact, d, vmultc, z};
        }
        cviol = std::max(0.0, -b.minCoeff());
        vmultc = Eigen::VectorXd::Constant(mcon, cviol) + b;
        z = Eigen::MatrixXd::Identity(num_vars, num_vars);

        // Check whether a quick return is possible
        if (cviol <= 0) {
            return {iact, nact, d, vmultc, z};
        }

        bool all_b_nan = true;
        for (int i = 0; i < b.size(); ++i)
            if (!std::isnan(b[i])) { all_b_nan = false; break; }
        if (all_b_nan) {
            return {iact, nact, d, vmultc, z};
        }

        // icon = index of most violated constraint (max(-b))
        icon = 0;
        double max_neg_b = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < mcon; ++i) {
            if (!std::isnan(b[i]) && -b[i] > max_neg_b) {
                max_neg_b = -b[i];
                icon = i;
            }
        }
        num_constraints = mcon;
    } else {
        // Stage 2: check whether a quick return is possible
        if (d.squaredNorm() >= delta * delta) {
            return {iact, nact, d, vmultc, z};
        }

        iact[mcon - 1] = mcon - 1;
        vmultc[mcon - 1] = 0;
        num_constraints = mcon - 1;
        icon = mcon - 1;

        // Recalculate cviol so it need not be passed from stage 1 to 2
        cviol = 0;
        for (int i = 0; i < num_constraints; ++i) {
            double adi = inprod(d, A.col(i));
            cviol = std::max(cviol, adi - b[i]);
        }
        cviol = std::max(0.0, cviol);
    }

    // zdota: scalar products of Z columns with active constraint gradients
    for (int k = 0; k < nact; ++k)
        zdota[k] = inprod(z.col(k), A.col(iact[k]));

    // More initialization
    double optold = REALMAX;
    int nactold = nact;
    int nfail = 0;
    int maxiter = std::min(10000, 100 * std::max(num_constraints, num_vars));

    // vmultd is computed from scratch at each iteration, but vmultc is inherited
    Eigen::VectorXd sdirn(num_vars);
    sdirn.setZero();

    for (int iter = 0; iter < maxiter; ++iter) {
        double optnew;
        if (stage == 1)
            optnew = cviol;
        else
            optnew = inprod(d, A.col(mcon - 1));

        // End the current stage if 3 consecutive iterations have failed to reduce
        // the best calculated objective value or to increase the number of active
        // constraints since the best value was calculated. This prevents cycling,
        // but there is a remote possibility of premature termination.
        if (optnew < optold || nact > nactold) {
            nactold = nact;
            nfail = 0;
        } else {
            nfail += 1;
        }
        optold = std::min(optold, optnew);
        if (nfail == 3) break;

        // If icon exceeds nact, add the constraint with index iact[icon] to the
        // active set.
        if (icon >= nact) {
            zdasav.head(nact) = zdota.head(nact);
            int nactsav = nact;

            Eigen::VectorXd new_col = A.col(iact[icon]);
            auto [z_new, zdota_new, nact_new] = qradd_Rdiag(new_col, z, zdota, nact);
            z = z_new;
            zdota = zdota_new;
            nact = nact_new;

            if (nact == nactsav + 1) {
                // qradd_Rdiag succeeded: the constraint was added to the active set
                if (nact != icon + 1) {
                    std::swap(vmultc[icon], vmultc[nact - 1]);
                    std::swap(iact[icon], iact[nact - 1]);
                } else {
                    vmultc[nact - 1] = 0;
                }
            } else {
                // qradd_Rdiag failed to add the constraint.
                // If the active set is empty, there is nothing to do.
                if (nact == 0) break;

                // VMULTD is calculated from scratch for the first (out of 2) time
                // in this iteration.
                // Note that IACT has not been updated to replace IACT[NACT] with
                // IACT[ICON]. Thus A[:, IACT[:NACT]] is the UNUPDATED version
                // before QRADD (note Z[:, :NACT] remains the same before and after
                // QRADD). Therefore if we supply ZDOTA to LSQR (as Rdiag) as Powell
                // did, we should use the UNUPDATED version, namely ZDASAV.
                Eigen::MatrixXd A_sub(num_vars, nact);
                for (int k = 0; k < nact; ++k) A_sub.col(k) = A.col(iact[k]);
                Eigen::VectorXd v = lsqr(A_sub, A.col(iact[icon]),
                                          z.leftCols(nact), zdasav.head(nact));
                vmultd.head(nact) = v;

                // N.B.: This can be triggered by NACT == 0 (among other
                // possibilities)! This is important, because NACT will be used as
                // an index in the sequel. However, the nact == 0 case is handled
                // above.
                bool any_positive = false;
                for (int k = 0; k < nact; ++k) {
                    if (vmultd[k] > 0 && iact[k] <= num_constraints) {
                        any_positive = true;
                        break;
                    }
                }
                if (!any_positive)
                    break;

                // vmultd[NACT+1:mcon] is not used, but initialize to avoid issues
                vmultd.tail(mcon - nact).setConstant(-1);

                // Revise the Lagrange multipliers.
                // Only the places with vmultd > 0 and iact <= num_constraints matter
                Eigen::VectorXd fracmult(nact);
                for (int k = 0; k < nact; ++k) {
                    if (vmultd[k] > 0 && iact[k] <= num_constraints)
                        fracmult[k] = vmultc[k] / vmultd[k];
                    else
                        fracmult[k] = REALMAX;
                }
                double frac = fracmult.minCoeff();

                for (int k = 0; k < nact; ++k)
                    vmultc[k] = std::max(0.0, vmultc[k] - frac * vmultd[k]);

                // Exit if the new value of zdota[nact] is not acceptable.
                // Powell's condition: not abs(zdota[nact]) > 0.
                // Note that it is different from 'abs(zdota[nact]) <= 0)' as
                // zdota[nact] can be NaN.
                // N.B.: We cannot arrive here with nact == 0, which should have
                // triggered a break above.
                if (std::isnan(zdota[nact - 1]) || std::abs(zdota[nact - 1]) <= EPS * EPS)
                    break;

                // Reorder the active constraints so that the one to be replaced is
                // at the end of the list.
                vmultc[icon] = 0;
                vmultc[nact - 1] = frac;
                std::swap(iact[icon], iact[nact - 1]);
            }

            // In stage 2, ensure that the objective continues to be treated as the
            // last active constraint.
            if (stage == 2 && iact[nact - 1] != (mcon - 1)) {
                if (nact <= 1) {
                    // We must exit, as nact-2 is used as an index below.
                    break;
                }
                Eigen::MatrixXd A_sub_full(num_vars, nact);
                for (int k = 0; k < nact; ++k) A_sub_full.col(k) = A.col(iact[k]);
                auto [z_exc, zdota_exc] = qrexc_Rdiag(A_sub_full, z, zdota.head(nact), nact - 2);
                z = z_exc;
                zdota.head(nact) = zdota_exc;
                std::swap(iact[nact - 2], iact[nact - 1]);
                std::swap(vmultc[nact - 2], vmultc[nact - 1]);
            }

            // Powell's code does not have the following. It avoids subsequent
            // floating point exceptions.
            if (std::isnan(zdota[nact - 1]) || std::abs(zdota[nact - 1]) <= EPS * EPS)
                break;

            // Set sdirn to the direction of the next change to the current vector.
            // During stage 1, sdirn gives a search direction that reduces all the
            // active constraint violations by one simultaneously.
            if (stage == 1) {
                sdirn -= ((inprod(sdirn, A.col(iact[nact - 1])) + 1) / zdota[nact - 1]) * z.col(nact - 1);
            } else {
                sdirn = (-1.0 / zdota[nact - 1]) * z.col(nact - 1);
            }
        } else {
            // icon < nact: delete the constraint with index iact[icon] from the
            // active set. Reorder iact[icon:nact] into [iact[icon+1:nact], iact[icon]]
            // then reduce nact to nact-1.
            Eigen::MatrixXd A_sub_full(num_vars, nact);
            for (int k = 0; k < nact; ++k) A_sub_full.col(k) = A.col(iact[k]);
            auto [z_exc, zdota_exc] = qrexc_Rdiag(A_sub_full, z, zdota.head(nact), icon);
            z = z_exc;
            zdota.head(nact) = zdota_exc;

            int saved_icon_val = iact[icon];
            int saved_vmultc_val = vmultc[icon];
            for (int k = icon; k < nact - 1; ++k) {
                iact[k] = iact[k + 1];
                vmultc[k] = vmultc[k + 1];
            }
            iact[nact - 1] = saved_icon_val;
            vmultc[nact - 1] = saved_vmultc_val;
            nact -= 1;

            // In theory, nact > 0 in stage 2, as the objective function should
            // always be considered as an "active constraint". However, looking at
            // the code, it's not guaranteed. It did happen in stage 1 that nact
            // became 0 after the reduction after almost one year of random tests.
            if (stage == 2 && nact < 0) break;
            // Powell's code does not have the following. Avoids exceptions.
            if (nact > 0) {
                if (std::isnan(zdota[nact - 1]) || std::abs(zdota[nact - 1]) <= EPS * EPS)
                    break;
            }

            // Set sdirn to the direction of the next change to the current vector.
            if (stage == 1) {
                if (nact < num_vars)
                    sdirn -= inprod(sdirn, z.col(nact)) * z.col(nact);
                // sdirn is orthogonal to z[:, nact+1]
            } else {
                sdirn = (-1.0 / zdota[nact - 1]) * z.col(nact - 1);
            }
        }

        // Calculate the step to the trust-region boundary or the step that reduces
        // cviol to 0.
        // The following calculation of step is adopted from NEWUOA/BOBYQA/LINCOA.
        // It seems to improve the performance of COBYLA. We also found that removing
        // the precaution about underflows is beneficial --- the underflows are
        // harmless anyway.
        double dd = delta * delta - d.squaredNorm();
        double ss = sdirn.squaredNorm();
        double sd = inprod(sdirn, d);
        if (dd <= 0 || ss <= EPS * delta * delta || std::isnan(sd)) break;

        // sqrtd: square root of a discriminant. The max avoids sqrtd < abs(sd) due
        // to underflow.
        double sqrtd = std::max(std::sqrt(ss * dd + sd * sd), std::abs(sd));
        sqrtd = std::max(sqrtd, std::sqrt(ss * dd));

        double step;
        if (sd > 0)
            step = dd / (sqrtd + sd);
        else
            step = (sqrtd - sd) / ss;

        // step < 0 should not happen. Step can be 0 or NaN when, e.g., sd or ss
        // becomes inf.
        if (step <= 0 || !std::isfinite(step)) break;

        if (stage == 1) {
            if (isminor(cviol, step)) break;
            step = std::min(step, cviol);
        }

        // Set dnew to the new variables if step is the steplength, and reduce cviol
        // to the corresponding maximum residual if stage 1 is being done.
        Eigen::VectorXd dnew = d + step * sdirn;

        if (stage == 1) {
            double max_res = 0;
            for (int k = 0; k < nact; ++k) {
                double res = inprod(dnew, A.col(iact[k])) - b[iact[k]];
                max_res = std::max(max_res, res);
            }
            cviol = std::max(0.0, max_res);
        }

        // vmultd is computed from scratch for the second (out of 2) time in one
        // iteration. vmultd[:nact] and vmultd[nact:mcon] are calculated separately
        // with no coupling. vmultd will be calculated from scratch again in the
        // next iteration.
        // Set vmultd to the vmultc vector that would occur if d became dnew.
        // A device is included to force vmultd[k] = 0 if deviations can be
        // attributed to computer rounding errors.
        Eigen::VectorXd vmultd_new(mcon);
        {
            Eigen::MatrixXd A_sub_nact(num_vars, nact);
            for (int k = 0; k < nact; ++k) A_sub_nact.col(k) = A.col(iact[k]);
            Eigen::VectorXd v = lsqr(A_sub_nact, dnew, z.leftCols(nact), zdota.head(nact));
            for (int k = 0; k < nact; ++k) vmultd_new[k] = -v[k];
        }
        if (stage == 2)
            vmultd_new[nact - 1] = std::max(0.0, vmultd_new[nact - 1]);

        // Complete vmultd by finding the new constraint residuals.
        Eigen::VectorXd cvshift(mcon);
        for (int k = 0; k < mcon; ++k)
            cvshift[k] = cviol - (inprod(dnew, A.col(iact[k])) - b[iact[k]]);
        Eigen::VectorXd cvsabs(mcon);
        for (int k = 0; k < mcon; ++k)
            cvsabs[k] = inprod(dnew.cwiseAbs(), A.col(iact[k]).cwiseAbs()) + std::abs(b[iact[k]]) + cviol;
        apply_isminor(cvshift, cvsabs);
        for (int k = nact; k < mcon; ++k)
            vmultd_new[k] = cvshift[k];

        // Calculate the fraction of the step from d to dnew that will be taken.
        double frac = 1.0;
        icon = -1;
        for (int k = 0; k < mcon; ++k) {
            if (vmultd_new[k] < 0) {
                double fk = vmultc[k] / (vmultc[k] - vmultd_new[k]);
                if (fk < frac) {
                    frac = fk;
                    icon = k;
                }
            }
        }

        // Update d, vmultc, and cviol
        Eigen::VectorXd dold = d;
        d = (1.0 - frac) * d + frac * dnew;
        for (int k = 0; k < mcon; ++k)
            vmultc[k] = std::max(0.0, (1.0 - frac) * vmultc[k] + frac * vmultd_new[k]);

        // Break in the case of inf/nan in d or vmultc
        if (!std::isfinite(d.cwiseAbs().sum()) || !std::isfinite(vmultc.cwiseAbs().sum())) {
            d = dold;
            break;
        }

        if (stage == 1) {
            double max_res2 = 0;
            for (int k = 0; k < mcon; ++k) {
                double res = inprod(d, A.col(k)) - b[k];
                max_res2 = std::max(max_res2, res);
            }
            // In theory, cviol = max(d@A - b, 0), yet the cviol updated as above
            // can be quite different from this value if A has huge entries (> 1e20).
            cviol = std::max(0.0, max_res2);
        }

        // In Powell's code, the condition is icon == 0. Indeed, icon < 0 cannot
        // hold unless fracmult contains only nan, which should not happen; icon >=
        // mcon should never occur.
        if (icon < 0 || icon >= mcon) break;
    }

    return {iact, nact, d, vmultc, z};
}

} // namespace detail

// TRSTLP: calculate an n-component vector d by the following two stages.
// Stage 1: d is set to the shortest vector that minimizes the greatest violation
// of the constraints A^T D <= B, subject to the Euclidean length of d being at
// most delta. If its length is strictly less than delta, then stage 2 uses the
// resultant freedom in d to minimize G^T D subject to no increase in any greatest
// constraint violation.
//
// cviol is the largest constraint violation of the current d: max(max(A^T D - b), 0).
// icon is the index of a most violated constraint if cviol is positive.
//
// nact is the number of constraints in the active set and iact[0..nact-1] are
// their indices, while the remainder of iact contains a permutation of the
// remaining constraint indices.
// N.B.: nact <= min(num_constraints, num_vars).
//
// Z is an orthogonal matrix whose first nact columns can be regarded as the
// result of Gram-Schmidt applied to the active constraint gradients. For
// j = 0..nact-1, zdota[j] is the scalar product of the jth column of Z with
// the gradient of the jth active constraint. d is the current vector of variables
// and here the residuals of the active constraints should be zero.
//
// N.B.:
// 0. In Powell's implementation, the constraints are A^T D >= B. In other words,
//    the A and B in our implementation are the negative of those in Powell's.
// 1. The algorithm was NOT documented in the COBYLA paper.
// 2. As a major part of the algorithm (see trstlp_sub), the code maintains and
//    updates the QR factorization of A[iact[:nact]], i.e. the gradients of all
//    active constraints. Z is indeed Q, and zdota is the diagonal of R.
// 3. There are probably better algorithms available for this problem.
inline Eigen::VectorXd trstlp(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                               double delta, const Eigen::VectorXd& g) {
    int num_constraints = A.cols();
    int num_vars = A.rows();

    Eigen::VectorXd vmultc;
    Eigen::VectorXi iact;
    int nact = 0;
    Eigen::VectorXd d = Eigen::VectorXd::Zero(num_vars);
    Eigen::MatrixXd z = Eigen::MatrixXd::Identity(num_vars, num_vars);

    if (num_constraints <= 0) {
        // No constraints: steepest descent direction, scaled to trust region
        double gn = norm(g);
        if (gn > 0) d = -(delta / gn) * g;
        return d;
    }

    vmultc = Eigen::VectorXd::Zero(num_constraints + 1);
    iact = Eigen::VectorXi::Zero(num_constraints + 1);

    // Form A_aug and b_aug. This allows the gradient of the objective to be
    // regarded as the gradient of a constraint in the second stage.
    Eigen::MatrixXd A_aug(num_vars, num_constraints + 1);
    A_aug.leftCols(num_constraints) = A;
    A_aug.rightCols(1).col(0) = g;
    Eigen::VectorXd b_aug(num_constraints + 1);
    b_aug.head(num_constraints) = b;
    b_aug[num_constraints] = 0;

    // Scale the problem if A contains large values. Otherwise floating point
    // exceptions may occur. Note that the trust-region step is scale invariant.
    for (int i = 0; i < num_constraints + 1; ++i) {
        double maxval = A_aug.col(i).cwiseAbs().maxCoeff();
        if (maxval > 1e12) {
            double modscal = std::max(2.0 * REALMIN, 1.0 / maxval);
            A_aug.col(i) *= modscal;
            b_aug[i] *= modscal;
        }
    }

    // Stage 1: minimize the L-infinity constraint violation of the linearized constraints
    {
        Eigen::VectorXi iact_sub = iact.head(num_constraints);
        Eigen::VectorXd vmultc_sub = vmultc.head(num_constraints);
        auto [iact_out, nact_out, d_out, vmultc_out, z_out] =
            detail::trstlp_sub(iact_sub, nact, 1,
                               A_aug.leftCols(num_constraints),
                               b_aug.head(num_constraints),
                               delta, d, vmultc_sub, z);
        iact.head(num_constraints) = iact_out;
        nact = nact_out;
        d = d_out;
        vmultc.head(num_constraints) = vmultc_out;
        z = z_out;
    }

    // Stage 2: minimize the linearized objective without increasing the L-infinity
    // constraint violation
    {
        auto [iact_out, nact_out, d_out, vmultc_out, z_out] =
            detail::trstlp_sub(iact, nact, 2, A_aug, b_aug, delta, d, vmultc, z);
        iact = iact_out;
        nact = nact_out;
        d = d_out;
        vmultc = vmultc_out;
        z = z_out;
    }

    return d;
}

// TRRAD: update the trust-region radius according to RATIO and DNORM.
inline double trrad(double delta_in, double dnorm, double eta1, double eta2,
                    double gamma1, double gamma2, double ratio) {
    double delta;
    if (ratio <= eta1) {
        // Powell's UOBYQA/NEWUOA
        // Powell's COBYLA/LINCOA: delta = gamma1 * delta_in
        // Powell's BOBYQA: delta = min(gamma1 * delta_in, dnorm)
        delta = gamma1 * dnorm;
    } else if (ratio <= eta2) {
        // Powell's UOBYQA/NEWUOA/BOBYQA/LINCOA
        delta = std::max(gamma1 * delta_in, dnorm);
    } else {
        // Powell's NEWUOA/BOBYQA
        // Modified version: delta = max(delta_in, gamma2 * dnorm). Works well for
        // UOBYQA. For noise-free CUTEst problems of <= 100 variables, Powell's
        // version works slightly better.
        // Powell's LINCOA: delta = min(max(gamma1*delta_in, gamma2*dnorm), gamma3*delta_in)
        delta = std::max(gamma1 * delta_in, gamma2 * dnorm);
    }
    return delta;
}

} // namespace prima

#endif // PRIMA_CPP_COBYLA_TRUSTREGION_HPP
