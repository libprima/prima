#ifndef PRIMA_CPP_COBYLA_COBYLB_HPP
#define PRIMA_CPP_COBYLA_COBYLB_HPP

// This module performs the major calculations of COBYLA.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
//
// N.B. (Zaikun 20220131): Powell's implementation of COBYLA uses RHO rather than
// DELTA as the trust-region radius, and RHO is never increased. DELTA does not exist
// in Powell's COBYLA code. Following the idea in Powell's other solvers (UOBYQA, ...,
// LINCOA), our code uses DELTA as the trust-region radius, while RHO works as a lower
// bound of DELTA and indicates the current resolution of the algorithm. DELTA is
// updated in a classical way subject to DELTA >= RHO, whereas RHO is updated as in
// Powell's COBYLA code and is never increased. The new implementation improves the
// performance of COBYLA.

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <functional>
#include <tuple>
#include <vector>

#include "../checkbreak.hpp"
#include "../consts.hpp"
#include "../evaluate.hpp"
#include "../history.hpp"
#include "../infos.hpp"
#include "../linalg.hpp"
#include "../message.hpp"
#include "../ratio.hpp"
#include "../redrho.hpp"
#include "../selectx.hpp"
#include "geometry.hpp"
#include "initialize.hpp"
#include "trustregion.hpp"
#include "update.hpp"

namespace prima {

// FC RATIO: calculate the ratio between the "typical change" of F and that of CONSTR.
// See equations (12)-(13) in Section 3 of the COBYLA paper for the definition.
inline double fcratio(const Eigen::MatrixXd& conmat, const Eigen::VectorXd& fval) {
    int m = conmat.rows();
    if (m == 0) return 0;
    Eigen::VectorXd cmin = -conmat.rowwise().maxCoeff();
    Eigen::VectorXd cmax = -conmat.rowwise().minCoeff();
    double fmin = fval.minCoeff();
    double fmax = fval.maxCoeff();
    double r = 0;
    if ((cmin.array() < 0.5 * cmax.array()).any() && fmin < fmax) {
        Eigen::VectorXd mask = (cmin.array() < 0.5 * cmax.array()).select(Eigen::VectorXd::Ones(m), Eigen::VectorXd::Zero(m));
        // Powell mentioned the following alternative in section 4 of his COBYLA paper.
        // According to a test on 20230610, it does not make much difference.
        // denom = np.max(max(cmax, 0) - cmin, mask=(cmin < 0.5 * cmax))
        double denom = std::numeric_limits<double>::infinity();
        for (int i = 0; i < m; ++i) {
            if (mask[i] > 0) {
                denom = std::min(denom, std::max(cmax[i], 0.0) - cmin[i]);
            }
        }
        if (std::isfinite(denom) && denom > 0)
            r = (fmax - fmin) / denom;
    }
    return r;
}

// GET CPEN: get the penalty parameter CPEN so that PREREM = PREREF + CPEN * PREREC > 0.
// See the discussions around equation (9) of the COBYLA paper.
//
// Increase CPEN if necessary to ensure PREREM > 0. Branch back for the next loop
// if this change alters the optimal vertex of the current simplex.
//
// Note the following:
// 1. In each loop, CPEN is changed only if PREREC > 0 > PREREF, in which case
//    PREREM is guaranteed positive after the update. Note that PREREC >= 0 and
//    max(PREREC, PREREF) > 0 in theory. If this holds numerically as well then CPEN
//    is not changed only if PREREC = 0 or PREREF >= 0, in which case PREREM is
//    currently positive, explaining why CPEN needs no update.
// 2. Even without an upper bound for the loop counter, the loop can occur at most
//    NUM_VARS+1 times. This is because the update of CPEN does not decrease CPEN,
//    and hence it can make vertex J (J <= NUM_VARS) become the new optimal vertex
//    only if CVAL[J] is less than CVAL[NUM_VARS], which can happen at most NUM_VARS
//    times. See the paragraph below (9) in the COBYLA paper. After the "correct"
//    optimal vertex is found, one more loop is needed to calculate CPEN, and hence
//    the loop can occur at most NUM_VARS+1 times.
// Powell's code defines BARMU = -PREREF / PREREC, and CPEN is increased to
// 2*BARMU if and only if it is currently less than 1.5*BARMU, a very
// "Powellful" scheme. In our implementation, however, we set CPEN directly to
// the maximum between its current value and 2*BARMU while handling possible
// overflow. This simplifies the scheme without worsening the performance.
inline double getcpen(const Eigen::MatrixXd* amat, const Eigen::VectorXd* bvec,
                      Eigen::MatrixXd& conmat, double cpen,
                      Eigen::VectorXd& cval, double delta,
                      Eigen::VectorXd& fval, double rho,
                      Eigen::MatrixXd& sim, Eigen::MatrixXd& simi) {
    int m_lcon = (bvec != nullptr) ? bvec->size() : 0;
    int num_constraints = conmat.rows();
    int num_vars = sim.rows();
    int info = INFO_DEFAULT;

    // Work on local copies to avoid modifying caller state (see pyprima bug HS102)
    Eigen::MatrixXd conmat_loc = conmat;
    Eigen::VectorXd cval_loc = cval;
    Eigen::VectorXd fval_loc = fval;
    Eigen::MatrixXd sim_loc = sim;
    Eigen::MatrixXd simi_loc = simi;

    Eigen::MatrixXd A(num_vars, num_constraints);

    for (int iter = 0; iter <= num_vars; ++iter) {
        // Switch the best vertex of the current simplex to SIM[:, NUM_VARS]
        std::tie(conmat_loc, cval_loc, fval_loc, sim_loc, simi_loc, info) =
            updatepole(cpen, conmat_loc, cval_loc, fval_loc, sim_loc, simi_loc);
        if (info == DAMAGING_ROUNDING) break;

        // Calculate the linear approximations to the objective and constraint functions
        Eigen::VectorXd fdiff1 = (fval_loc.head(num_vars).array() - fval_loc[num_vars]).matrix();
        Eigen::VectorXd g = matprod(fdiff1, simi_loc);
        A.setZero();
        if (amat != nullptr && m_lcon > 0) {
            Eigen::MatrixXd amat_t = amat->transpose();
            A.leftCols(m_lcon) = amat_t;
        }
        int m_nlcon = num_constraints - m_lcon;
        if (m_nlcon > 0) {
            Eigen::VectorXd col_nlcon = conmat_loc.col(num_vars).segment(m_lcon, m_nlcon);
            Eigen::MatrixXd diff = conmat_loc.block(m_lcon, 0, m_nlcon, num_vars).colwise()
                                   - col_nlcon;
            A.rightCols(m_nlcon) = matprod(diff, simi_loc).transpose();
        }

        // Calculate the trust-region trial step D. Note that D does NOT depend on CPEN.
        Eigen::VectorXd d = trstlp(A, -conmat_loc.col(num_vars), delta, g);

        // Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D
        double preref = -inprod(d, g);  // Can be negative
        double ad_max = 0;
        for (int i = 0; i < num_constraints; ++i)
            ad_max = std::max(ad_max, conmat_loc(i, num_vars) + inprod(d, A.col(i)));
        double prerec = cval_loc[num_vars] - std::max(0.0, ad_max);

        // PREREC <= 0 or PREREF >= 0 or either is NaN
        if (!(prerec > 0 && preref < 0)) break;

        cpen = std::max(cpen, std::min(-2 * preref / prerec, REALMAX));

        if (findpole(cpen, cval_loc, fval_loc) == num_vars) break;
    }

    return cpen;
}

inline std::tuple<Eigen::VectorXd, double, Eigen::VectorXd, double, int,
                  std::vector<Eigen::VectorXd>, std::vector<double>,
                  std::vector<double>, std::vector<Eigen::VectorXd>, int>
cobylb(const Calcfc& calcfc, int iprint, int maxfilt, int maxfun,
       const Eigen::MatrixXd* amat, const Eigen::VectorXd* bvec,
       double ctol, double cweight, double eta1, double eta2,
       double ftarget, double gamma1, double gamma2,
       double rhobeg, double rhoend,
       const Eigen::VectorXd& constr_in, double f_in,
       const Eigen::VectorXd& x_in, int maxhist,
       const std::function<bool(const Eigen::VectorXd&, double, int, int, double,
                                 const Eigen::VectorXd&)>& callback) {

    // Outputs
    std::vector<Eigen::VectorXd> xhist, conhist;
    std::vector<double> fhist, chist;

    // Sizes
    int m_lcon = (bvec != nullptr) ? bvec->size() : 0;
    int num_constraints = constr_in.size();
    int m_nlcon = num_constraints - m_lcon;
    int num_vars = x_in.size();

    // Initialize SIM, FVAL, CONMAT, and CVAL, together with the history.
    // After the initialization, SIM[:, NUM_VARS] holds the vertex of the initial
    // simplex with the smallest function value (regardless of the constraint
    // violation), and SIM[:, :NUM_VARS] holds the displacements from the other
    // vertices to SIM[:, NUM_VARS]. FVAL, CONMAT, and CVAL hold the function values,
    // constraint values, and constraint violations on the vertices in the order
    // corresponding to SIM.
    auto [evaluated, conmat, cval, sim, simi, fval, nf, subinfo] =
        initxfc(calcfc, iprint, maxfun, constr_in, amat, bvec, ctol, f_in,
                ftarget, rhobeg, x_in, xhist, fhist, chist, conhist, maxhist);

    // Initialize the filter, including xfilt, ffilt, confilt, cfilt, and nfilt.
    // N.B.: The filter is used only when selecting which iterate to return. It does
    // not interfere with the iterations. COBYLA is NOT a filter method but a
    // trust-region method based on an L-infinity merit function. Powell's
    // implementation does not use a filter to select the iterate, possibly returning
    // a suboptimal iterate.
    int maxfilt_actual = std::min(std::max(maxfilt, 1), maxfun);
    Eigen::VectorXd cfilt = Eigen::VectorXd::Zero(maxfilt_actual);
    Eigen::MatrixXd confilt = Eigen::MatrixXd::Zero(num_constraints, maxfilt_actual);
    Eigen::VectorXd ffilt = Eigen::VectorXd::Zero(maxfilt_actual);
    Eigen::MatrixXd xfilt = Eigen::MatrixXd::Zero(num_vars, maxfilt_actual);

    int nfilt = initfilt(conmat, ctol, cweight, cval, fval, sim, evaluated,
                         cfilt, confilt, ffilt, xfilt);

    // Check whether to return due to abnormal cases that may occur during init
    if (subinfo != INFO_DEFAULT) {
        // N.B.: selectx and findpole choose X by different standards, one cannot
        // replace the other
        int kopt = selectx(ffilt.head(nfilt), cfilt.head(nfilt), cweight, ctol);
        Eigen::VectorXd x = xfilt.col(kopt);
        double f = ffilt[kopt];
        Eigen::VectorXd constr = confilt.col(kopt);
        double cstrv = cfilt[kopt];
        retmsg("COBYLA", subinfo, iprint, nf, f, x, cstrv, &constr);
        return {x, f, constr, cstrv, nf, xhist, fhist, chist, conhist, subinfo};
    }

    // Set some more initial values.
    // We must initialize shortd, ratio, and jdrop_tr because these get defined on
    // branches that are not guaranteed to be executed, but their values are used later.
    // Our initialization of CPEN differs from Powell's in two ways. First, we use the
    // ratio defined in (13) of Powell's COBYLA paper to initialize CPEN. Second, we
    // impose CPEN >= CPENMIN > 0. Powell's code simply initializes CPEN to 0.
    double rho = rhobeg;
    double delta = rhobeg;
    double cpenmin = EPS;
    double cpen = std::max(cpenmin, std::min(1.0e3, fcratio(conmat, fval)));
    bool shortd = false;
    double ratio = -1;
    int jdrop_tr = 0;
    Eigen::VectorXd d_last = Eigen::VectorXd::Zero(num_vars);

    // If DELTA <= GAMMA3*RHO after an update, we set DELTA to RHO. GAMMA3 must be
    // less than GAMMA2. The reason:
    // Imagine a very successful step with DNORM = the un-updated DELTA = RHO. TRRAD
    // will update DELTA to GAMMA2*RHO. If GAMMA3 >= GAMMA2, then DELTA will be reset
    // to RHO, which is not reasonable as D is very successful. See paragraph two of
    // Sec 5.2.5 in T. M. Ragonneau's thesis: "Model-Based Derivative-Free Optimization
    // Methods and Software."
    double gamma3 = std::max(1.0, std::min(0.75 * gamma2, 1.5));

    // MAXTR is the maximal number of trust-region iterations. Each trust-region
    // iteration takes 1 or 2 function evaluations unless the trust-region step is
    // short or the trust-region subproblem solver fails but the geometry step is not
    // invoked. Thus the following MAXTR is unlikely to be reached.
    int maxtr = 10 * maxfun;
    int info = MAXTR_REACHED;

    Eigen::MatrixXd A(num_vars, num_constraints);
    Eigen::VectorXd distsq(num_vars + 1);

    // Begin the iterative procedure
    // After solving a trust-region subproblem, we use three boolean variables to
    // control the workflow.
    // SHORTD - Is the trust-region trial step too short to invoke a function eval?
    // IMPROVE_GEO - Will we improve the model after the trust-region iteration? If
    //               yes, a geometry step will be taken, corresponding to "Branch
    //               (Delta)" in the COBYLA paper.
    // REDUCE_RHO - Will we reduce rho after the trust-region iteration?
    // COBYLA never sets IMPROVE_GEO and REDUCE_RHO to True simultaneously.
    for (int tr = 0; tr < maxtr; ++tr) {
        // Increase the penalty parameter CPEN, if needed, so that
        // PREREM = PREREF + CPEN * PREREC > 0.
        // This is the first (out of two) update of CPEN, where CPEN increases or
        // remains the same.
        // N.B.: CPEN and the merit function PHI = FVAL + CPEN*CVAL are used in three
        // places only.
        // 1. In FINDPOLE/UPDATEPOLE, deciding the optimal vertex of the simplex.
        // 2. After the trust-region trial step, calculating the reduction ratio.
        // 3. In GEOSTEP, deciding the direction of the geometry step.
        // They do not appear explicitly in the trust-region subproblem, though the
        // trust-region center (i.e. the current optimal vertex) is defined by them.
        cpen = getcpen(amat, bvec, conmat, cpen, cval, delta, fval, rho, sim, simi);

        // Switch the best vertex of the current simplex to SIM[:, NUM_VARS]
        std::tie(conmat, cval, fval, sim, simi, subinfo) =
            updatepole(cpen, conmat, cval, fval, sim, simi);
        if (subinfo == DAMAGING_ROUNDING) { info = subinfo; break; }

        // Does the interpolation set have adequate geometry? It affects improve_geo
        // and reduce_rho.
        bool adequate_geo = true;
        for (int i = 0; i < num_vars; ++i) {
            double nsq = sim.col(i).squaredNorm();
            if (nsq > 4 * delta * delta) { adequate_geo = false; break; }
        }

        // Calculate the linear approximations to the objective and constraint functions.
        // N.B.: TRSTLP accesses A mostly by columns, so it is more reasonable to
        // store A instead of A^T.
        Eigen::VectorXd fdiff2 = (fval.head(num_vars).array() - fval[num_vars]).matrix();
        Eigen::VectorXd g = matprod(fdiff2, simi);
        A.setZero();
        if (amat != nullptr && m_lcon > 0) {
            Eigen::MatrixXd amat_t = amat->transpose();
            A.leftCols(m_lcon) = amat_t;
        }
        int mnlcon = num_constraints - m_lcon;
        if (mnlcon > 0) {
            Eigen::VectorXd col_nlcon2 = conmat.col(num_vars).segment(m_lcon, mnlcon);
            Eigen::MatrixXd diff = conmat.block(m_lcon, 0, mnlcon, num_vars).colwise()
                                   - col_nlcon2;
            A.rightCols(mnlcon) = matprod(diff, simi).transpose();
        }

        // Calculate the trust-region trial step d. Note that d does NOT depend on cpen.
        d_last = trstlp(A, -conmat.col(num_vars), delta, g);
        Eigen::VectorXd d = d_last;
        double dnorm = std::min(delta, norm(d));

        // Is the trust-region trial step short? N.B.: we compare DNORM with RHO, not
        // DELTA. Powell's code especially defines SHORTD by SHORTD = (DNORM < 0.5 *
        // RHO). In our tests 1/10 seems to work better than 1/2 or 1/4, especially
        // for linearly constrained problems.
        shortd = (dnorm <= 0.1 * rho);

        // Predict the change to F (PREREF) and to the constraint violation (PREREC)
        // due to D. In theory:
        // 1. B[:NUM_CONSTRAINTS] = -CONMAT[:, NUM_VARS] and hence
        //    max(B[:NUM_CONSTRAINTS] - D@A[:, :NUM_CONSTRAINTS], 0) is the L-infinity
        //    violation of the linearized constraints. When D=0, the violation is
        //    CVAL[NUM_VARS]. PREREC is the reduction of this violation achieved by D,
        //    which is nonnegative in theory; PREREC = 0 iff B[:NUM_CONSTRAINTS] <= 0,
        //    i.e. the trust-region center satisfies the linearized constraints.
        // 2. PREREF may be negative or 0, but it is positive when PREREC = 0 and
        //    shortd is False.
        // 3. Due to 2, in theory, max(PREREC, PREREF) > 0 if shortd is False.
        double preref = -inprod(d, g);  // Can be negative
        double ad_max = 0;
        for (int i = 0; i < num_constraints; ++i)
            ad_max = std::max(ad_max, conmat(i, num_vars) + inprod(d, A.col(i)));
        double prerec = cval[num_vars] - std::max(0.0, ad_max);

        // Evaluate PREREM, the predicted reduction in the merit function.
        // In theory, PREREM >= 0 and it is 0 iff CPEN = 0 = PREREF.
        double prerem = preref + cpen * prerec;
        bool trfail = !(prerem > 1.0e-6 * std::min(cpen, 1.0) * rho);

        if (shortd || trfail) {
            // Reduce DELTA if D is short or fails to render PREREM > 0. The latter
            // can only happen due to rounding errors.
            delta *= 0.1;
            if (delta <= gamma3 * rho) delta = rho;
        } else {
            // Calculate the next value of the objective and constraint functions.
            // If X is close to one of the points in the interpolation set, then we
            // do not evaluate, assuming the values at the closest point.
            // N.B.: If this happens, do NOT include X into the filter, as F and
            // CONSTR are inaccurate.
            Eigen::VectorXd x = sim.col(num_vars) + d;
            distsq[num_vars] = (x - sim.col(num_vars)).squaredNorm();
            for (int i = 0; i < num_vars; ++i)
                distsq[i] = (x - (sim.col(num_vars) + sim.col(i))).squaredNorm();

            int j;
            distsq.minCoeff(&j);
            double f; Eigen::VectorXd constr; double cstrv;
            double rhoend_sq = 1e-4 * rhoend;
            rhoend_sq = rhoend_sq * rhoend_sq;

            if (distsq[j] <= rhoend_sq) {
                f = fval[j]; constr = conmat.col(j); cstrv = cval[j];
            } else {
                std::tie(f, constr) = evaluate(calcfc, x, m_nlcon, amat, bvec);
                cstrv = 0;
                for (int i = 0; i < constr.size(); ++i) cstrv = std::max(cstrv, constr[i]);
                ++nf;
                savehist(maxhist, x, xhist, f, fhist, cstrv, chist, constr, conhist);
                nfilt = savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, &constr, &confilt);
            }

            fmsg("COBYLA", "Trust region", iprint, nf, delta, f, x, cstrv, &constr);

            // Evaluate ACTREM, the actual reduction in the merit function
            double actrem = (fval[num_vars] + cpen * cval[num_vars]) - (f + cpen * cstrv);

            // Calculate the reduction ratio by redrat, which handles inf/nan
            ratio = redrat(actrem, prerem, eta1);

            // Update DELTA. After this, DELTA < DNORM may hold.
            // N.B.:
            // 1. Powell's code uses RHO as the trust-region radius and updates it as:
            //    Reduce RHO to GAMMA1*RHO if ADEQUATE_GEO is TRUE and either SHORTD
            //    is TRUE or RATIO < ETA1, then revise RHO to RHOEND if its new value
            //    is not more than GAMMA3*RHOEND; RHO is never increased.
            // 2. Our implementation uses DELTA as the trust-region radius, with RHO
            //    as a lower bound for DELTA. DELTA is updated in a typical trust-
            //    region way, and revised to RHO if its new value is not more than
            //    GAMMA3*RHO. RHO reflects the current resolution of the algorithm.
            // 3. The same as Powell's code, we do not reduce RHO unless ADEQUATE_GEO
            //    is TRUE. This is also how Powell updated RHO in UOBYQA/NEWUOA/
            //    BOBYQA/LINCOA.
            delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio);
            if (delta <= gamma3 * rho) delta = rho;

            // Is the newly generated X better than the current best point?
            bool ximproved = actrem > 0;

            // Set JDROP_TR to the index of the vertex to be replaced with X.
            // JDROP_TR = -1 means there is no good point to replace, and X will not
            // be included into the simplex; in this case, the geometry likely needs
            // improvement, which will be handled below.
            jdrop_tr = setdrop_tr(ximproved, d, delta, rho, sim, simi);

            // Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM[:, JDROP_TR] is
            // replaced with D. UPDATEXFC does nothing if JDROP_TR is -1.
            std::tie(sim, simi, fval, conmat, cval, subinfo) =
                updatexfc(jdrop_tr, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi);
            if (subinfo == DAMAGING_ROUNDING) { info = subinfo; break; }

            // Check whether to break due to maxfun, ftarget, etc.
            subinfo = checkbreak_con(maxfun, nf, cstrv, ctol, f, ftarget, x);
            if (subinfo != INFO_DEFAULT) { info = subinfo; break; }
        }
        // End of if SHORTD or TRFAIL. The normal trust-region calculation ends.

        // Before the next trust-region iteration, we possibly improve the geometry
        // of the simplex or reduce RHO according to IMPROVE_GEO and REDUCE_RHO.
        // N.B.: We must ensure that the algorithm does not set IMPROVE_GEO = True
        // at infinitely many consecutive iterations without moving SIM[:, NUM_VARS]
        // or reducing RHO. This is ensured by:
        // 1. If an iteration sets IMPROVE_GEO to True, it must also reduce DELTA or
        //    set DELTA to RHO.
        // 2. If SIM[:, NUM_VARS] and RHO remain unchanged, then ADEQUATE_GEO will
        //    become True after at most NUM_VARS invocations of GEOSTEP.

        // BAD_TRSTEP: Is the last trust-region step bad?
        bool bad_trstep = shortd || trfail || ratio <= 0 || jdrop_tr == -1;
        // IMPROVE_GEO: Should we take a geometry step to improve geometry?
        bool improve_geo = bad_trstep && !adequate_geo;
        // REDUCE_RHO: Should we enhance the resolution by reducing rho?
        bool reduce_rho = bad_trstep && adequate_geo && std::max(delta, dnorm) <= rho;

        // COBYLA never sets IMPROVE_GEO and REDUCE_RHO to True simultaneously.

        // Improve the geometry of the simplex by removing a point and adding a new one.
        // If the current interpolation set has acceptable geometry, skip the step.
        if (improve_geo) {
            bool geo_ok = true;
            for (int i = 0; i < num_vars; ++i)
                if (sim.col(i).squaredNorm() > 4 * delta * delta) { geo_ok = false; break; }
            if (!geo_ok) {
                // Decide a vertex to drop from the simplex. It will be replaced with
                // SIM[:, NUM_VARS] + D to improve geometry.
                // N.B.:
                // 1. COBYLA never sets JDROP_GEO = num_vars.
                // 2. The following JDROP_GEO comes from UOBYQA/NEWUOA/BOBYQA/LINCOA.
                double max_dist = 0;
                int jdrop_geo = 0;
                for (int i = 0; i < num_vars; ++i) {
                    double nd = sim.col(i).squaredNorm();
                    if (nd > max_dist) { max_dist = nd; jdrop_geo = i; }
                }

                // Calculate the geometry step D
                double delbar = delta / 2;
                d = geostep(jdrop_geo, amat, bvec, conmat, cpen, cval, delbar, fval, simi);
                Eigen::VectorXd x = sim.col(num_vars) + d;

                distsq[num_vars] = (x - sim.col(num_vars)).squaredNorm();
                for (int i = 0; i < num_vars; ++i)
                    distsq[i] = (x - (sim.col(num_vars) + sim.col(i))).squaredNorm();
                int j; distsq.minCoeff(&j);
                double f; Eigen::VectorXd constr; double cstrv;
                double rhoend_sq = 1e-4 * rhoend;
                rhoend_sq = rhoend_sq * rhoend_sq;

                if (distsq[j] <= rhoend_sq) {
                    f = fval[j]; constr = conmat.col(j); cstrv = cval[j];
                } else {
                    std::tie(f, constr) = evaluate(calcfc, x, m_nlcon, amat, bvec);
                    cstrv = 0;
                    for (int i = 0; i < constr.size(); ++i) cstrv = std::max(cstrv, constr[i]);
                    ++nf;
                    savehist(maxhist, x, xhist, f, fhist, cstrv, chist, constr, conhist);
                    nfilt = savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, &constr, &confilt);
                }

                fmsg("COBYLA", "Geometry", iprint, nf, delta, f, x, cstrv, &constr);
                std::tie(sim, simi, fval, conmat, cval, subinfo) =
                    updatexfc(jdrop_geo, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi);
                if (subinfo == DAMAGING_ROUNDING) { info = subinfo; break; }

                subinfo = checkbreak_con(maxfun, nf, cstrv, ctol, f, ftarget, x);
                if (subinfo != INFO_DEFAULT) { info = subinfo; break; }
            }
        }

        // The calculations with the current RHO are complete. Enhance the resolution
        // by reducing RHO; update DELTA and CPEN at the same time.
        if (reduce_rho) {
            if (rho <= rhoend) { info = SMALL_TR_RADIUS; break; }
            delta = std::max(0.5 * rho, redrho(rho, rhoend));
            rho = redrho(rho, rhoend);
            // The second (out of two) update of CPEN, where CPEN decreases or
            // remains the same.
            // Powell's code: cpen = min(cpen, fcratio(fval, conmat)), may set to 0.
            cpen = std::max(cpenmin, std::min(cpen, fcratio(conmat, fval)));
            Eigen::VectorXd pole_constr = conmat.col(num_vars);
            rhomsg("COBYLA", iprint, nf, fval[num_vars], rho, sim.col(num_vars),
                   cval[num_vars], &pole_constr, cpen);
            std::tie(conmat, cval, fval, sim, simi, subinfo) =
                updatepole(cpen, conmat, cval, fval, sim, simi);
            if (subinfo == DAMAGING_ROUNDING) { info = subinfo; break; }
        }

        // Report the current best value, and check if user asks for early termination
        if (callback) {
            bool terminate = callback(sim.col(num_vars), fval[num_vars], nf, tr,
                                      cval[num_vars], conmat.col(num_vars));
            if (terminate) { info = CALLBACK_TERMINATE; break; }
        }
    }
    // End of for loop. The iterative procedure ends.

    // Return from the calculation, after trying the last trust-region step if it
    // has not been tried yet.
    // Ensure that D has not been updated after SHORTD == TRUE occurred, or the
    // code below is incorrect.
    Eigen::VectorXd x_try = sim.col(num_vars) + d_last;
    if (info == SMALL_TR_RADIUS && shortd &&
        (x_try - sim.col(num_vars)).norm() > 1.0e-3 * rhoend && nf < maxfun) {
        // UPDATEXFC or UPDATEPOLE is not called since the last trust-region step.
        // Hence SIM[:, NUM_VARS] remains unchanged. Otherwise SIM[:, NUM_VARS] + D
        // would not make sense.
        auto [f_try, constr_try] = evaluate(calcfc, x_try, m_nlcon, amat, bvec);
        double cstrv_try = 0;
        for (int i = 0; i < constr_try.size(); ++i) cstrv_try = std::max(cstrv_try, constr_try[i]);
        ++nf;
        savehist(maxhist, x_try, xhist, f_try, fhist, cstrv_try, chist, constr_try, conhist);
        nfilt = savefilt(cstrv_try, ctol, cweight, f_try, x_try, nfilt, cfilt, ffilt, xfilt, &constr_try, &confilt);
        // DELTA has been updated. RHO is only indicative here.
        fmsg("COBYLA", "Trust region", iprint, nf, rho, f_try, x_try, cstrv_try, &constr_try);
    }

    // Return the best calculated values of the variables
    // N.B.: SELECTX and FINDPOLE choose X by different standards, one cannot replace
    // the other.
    int kopt = selectx(ffilt.head(nfilt), cfilt.head(nfilt), std::max(cpen, cweight), ctol);
    Eigen::VectorXd x_ret = xfilt.col(kopt);
    double f_ret = ffilt[kopt];
    Eigen::VectorXd constr_ret = confilt.col(kopt);
    double cstrv_ret = cfilt[kopt];

    retmsg("COBYLA", info, iprint, nf, f_ret, x_ret, cstrv_ret, &constr_ret);

    return {x_ret, f_ret, constr_ret, cstrv_ret, nf,
            xhist, fhist, chist, conhist, info};
}

} // namespace prima

#endif // PRIMA_CPP_COBYLA_COBYLB_HPP
