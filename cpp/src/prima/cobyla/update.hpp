//----------------------------------------------------------------------------------------------------!
// This module contains subroutines concerning the update of the interpolation set.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
//
// Python translation by Nickolai Belakovski.
// C++ translation by ...
//----------------------------------------------------------------------------------------------------!

#ifndef PRIMA_CPP_COBYLA_UPDATE_HPP
#define PRIMA_CPP_COBYLA_UPDATE_HPP

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <tuple>

#include "../consts.hpp"
#include "../infos.hpp"
#include "../linalg.hpp"

namespace prima {

inline int findpole(double cpen, const Eigen::VectorXd& cval, const Eigen::VectorXd& fval) {
    // This subroutine identifies the best vertex of the current simplex with respect to the merit
    // function PHI = F + CPEN * CSTRV.
    int num_vars = fval.size() - 1;
    int jopt = num_vars;
    Eigen::VectorXd phi = fval + cpen * cval;
    double phimin = phi.minCoeff();

    // Identify the optimal vertex of the current simplex
    // Essentially jopt = argmin(phi). However, we keep jopt = num_vars unless there
    // is a strictly better choice. When there are multiple choices, we choose the jopt
    // with the smallest value of cval.
    if (phimin < phi[jopt] || ((cval.array() < cval[jopt]).any() && (phi.array() <= phi[jopt]).any())) {
        // While we could use argmin(phi), there may be two places where phi achieves
        // phimin, and in that case we should choose the one with the smallest cval.
        double best_cv = std::numeric_limits<double>::infinity();
        for (int i = 0; i <= num_vars; ++i) {
            if (phi[i] <= phimin && cval[i] < best_cv) {
                best_cv = cval[i];
                jopt = i;
            }
        }
    }
    return jopt;
}

inline std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd,
                  Eigen::MatrixXd, Eigen::MatrixXd, int>
updatepole(double cpen, const Eigen::MatrixXd& conmat, const Eigen::VectorXd& cval,
           const Eigen::VectorXd& fval, const Eigen::MatrixXd& sim,
           const Eigen::MatrixXd& simi) {
    //----------------------------------------------------------------------------------------------!
    // This subroutine identifies the best vertex of the current simplex with respect to the merit
    // function PHI = F + CPEN * CSTRV, and then switch this vertex to SIM[:, NUM_VARS], which Powell
    // called the "pole position" in his comments. CONMAT, CVAL, FVAL, and SIMI are updated
    // accordingly.
    //
    // N.B. 1: In precise arithmetic, the following two procedures produce the same results:
    // 1) apply UPDATEPOLE to SIM twice, first with CPEN = CPEN1 and then with CPEN = CPEN2;
    // 2) apply UPDATEPOLE to SIM with CPEN = CPEN2.
    // In finite-precision arithmetic, however, they may produce different results unless
    // CPEN1 = CPEN2.
    //
    // N.B. 2: When JOPT == N+1, the best vertex is already at the pole position, so there is nothing
    // to switch. However, as in Powell's code, the code below will check whether SIMI is good enough
    // to work as the inverse of SIM(:, 1:N) or not. If not, Powell's code would invoke an error
    // return of COBYLB; our implementation, however, will try calculating SIMI from scratch; if the
    // recalculated SIMI is still of poor quality, then UPDATEPOLE will return with INFO =
    // DAMAGING_ROUNDING, informing COBYLB that SIMI is poor due to damaging rounding errors.
    //
    // N.B. 3: UPDATEPOLE should be called when and only when FINDPOLE can potentially returns a value
    // other than N+1. The value of FINDPOLE is determined by CPEN, CVAL, and FVAL, the latter two
    // being decided by SIM. Thus UPDATEPOLE should be called after CPEN or SIM changes. COBYLA
    // updates CPEN at only two places: the beginning of each trust-region iteration, and when REDRHO
    // is called; SIM is updated only by UPDATEXFC, which itself calls UPDATEPOLE internally.
    // Therefore, we only need to call UPDATEPOLE after updating CPEN at the beginning of each
    // trust-region iteration and after each invocation of REDRHO.
    //----------------------------------------------------------------------------------------------!

    int num_constraints = conmat.rows();
    int num_vars = sim.rows();
    // INFO must be set, as it is an output.
    int info = INFO_DEFAULT;

    // Identify the optimal vertex of the current simplex.
    int jopt = findpole(cpen, cval, fval);

    Eigen::MatrixXd sim_new = sim;
    Eigen::MatrixXd simi_new = simi;
    Eigen::VectorXd fval_new = fval;
    Eigen::MatrixXd conmat_new = conmat;
    Eigen::VectorXd cval_new = cval;

    // Switch the best vertex to the pole position SIM[:, NUM_VARS] if it is not there already and
    // update SIMI. Before the update, save a copy of SIM and SIMI. If the update is unsuccessful due
    // to damaging rounding errors, we restore them and return with INFO = DAMAGING_ROUNDING.
    Eigen::MatrixXd sim_old = sim;
    Eigen::MatrixXd simi_old = simi;

    if (jopt >= 0 && jopt < num_vars) {
        // Unless there is a bug in FINDPOLE it is guaranteed that JOPT >= 0
        // When JOPT == NUM_VARS, there is nothing to switch; in addition SIMI[JOPT, :] will be
        // illegal.
        sim_new.col(num_vars) += sim_new.col(jopt);
        Eigen::VectorXd sim_jopt = sim_new.col(jopt);
        sim_new.col(jopt).setZero();
        for (int i = 0; i < num_vars; ++i)
            sim_new.col(i) -= sim_jopt;

        // The above update is equivalent to multiplying SIM[:, :NUM_VARS] from the right side by a
        // matrix whose JOPT-th row is [-1, -1, ..., -1], while all the other rows are the same as
        // those of the identity matrix. It is easy to check that the inverse of this matrix is
        // itself. Therefore, SIMI should be updated by a multiplication with this matrix (i.e. its
        // inverse) from the left side, as is done in the following line. The JOPT-th row of the
        // updated SIMI is minus the sum of all rows of the original SIMI, whereas all the other rows
        // remain unchanged.
        simi_new.row(jopt) = -simi.colwise().sum();
    }

    // Check whether SIMI is a poor approximation to the inverse of SIM[:, :NUM_VARS]
    // Calculate SIMI from scratch if the current one is damaged by rounding errors.
    Eigen::MatrixXd sim_left = sim_new.leftCols(num_vars);
    double erri = (matprod(simi_new, sim_left) -
                   Eigen::MatrixXd::Identity(num_vars, num_vars)).cwiseAbs().maxCoeff();
    double itol = 1;
    if (erri > 0.1 * itol || std::isnan(erri)) {
        Eigen::MatrixXd simi_test = inv(sim_left);
        double erri_test = (matprod(simi_test, sim_left) -
                            Eigen::MatrixXd::Identity(num_vars, num_vars)).cwiseAbs().maxCoeff();
        if (erri_test < erri || (std::isnan(erri) && !std::isnan(erri_test))) {
            simi_new = simi_test;
            erri = erri_test;
        }
    }

    // If SIMI is satisfactory, then update FVAL, CONMAT, and CVAL. Otherwise restore SIM and SIMI,
    // and return with INFO = DAMAGING_ROUNDING.
    if (erri <= itol) {
        if (jopt >= 0 && jopt < num_vars) {
            std::swap(fval_new[jopt], fval_new[num_vars]);
            conmat_new.col(jopt).swap(conmat_new.col(num_vars));
            std::swap(cval_new[jopt], cval_new[num_vars]);
        }
    } else {
        info = DAMAGING_ROUNDING;
        sim_new = sim_old;
        simi_new = simi_old;
    }

    return {conmat_new, cval_new, fval_new, sim_new, simi_new, info};
}

inline std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd,
                  Eigen::MatrixXd, Eigen::VectorXd, int>
updatexfc(int jdrop, const Eigen::VectorXd& constr, double cpen, double cstrv,
          const Eigen::VectorXd& d, double f, const Eigen::MatrixXd& conmat,
          const Eigen::VectorXd& cval, const Eigen::VectorXd& fval,
          const Eigen::MatrixXd& sim, const Eigen::MatrixXd& simi) {
    // This function revises the simplex by updating the elements of SIM, SIMI, FVAL, CONMAT, and
    // CVAL.

    int num_constraints = conmat.rows();
    int num_vars = sim.rows();
    int info = INFO_DEFAULT;

    // Do nothing when JDROP < 0. This can only happen after a trust-region step.
    // JDROP < 0 is impossible if the input is correct.
    if (jdrop < 0) {
        return {sim, simi, fval, conmat, cval, info};
    }

    Eigen::MatrixXd sim_new = sim;
    Eigen::MatrixXd simi_new = simi;

    if (jdrop < num_vars) {
        sim_new.col(jdrop) = d;
        double inprod_val = inprod(simi_new.row(jdrop), d);
        Eigen::VectorXd simi_jdrop = simi_new.row(jdrop).transpose() / inprod_val;
        simi_new -= outprod(matprod(simi_new, d), simi_jdrop);
        simi_new.row(jdrop) = simi_jdrop.transpose();
    } else {
        sim_new.col(num_vars) += d;
        for (int i = 0; i < num_vars; ++i)
            sim_new.col(i) -= d;

        Eigen::VectorXd simid = matprod(simi_new, d);
        double sum_simid = simid.sum();
        Eigen::VectorXd sum_simi = primasum(simi_new, 0);
        simi_new += outprod(simid, sum_simi / (1.0 - sum_simid));
    }

    // Check whether SIMI is a poor approximation to the inverse of SIM[:, :NUM_VARS]
    // Calculate SIMI from scratch if the current one is damaged by rounding errors.
    double itol = 1;
    Eigen::MatrixXd left = sim_new.leftCols(num_vars);
    double erri = (matprod(simi_new, left) -
                   Eigen::MatrixXd::Identity(num_vars, num_vars)).cwiseAbs().maxCoeff();
    if (erri > 0.1 * itol || std::isnan(erri)) {
        Eigen::MatrixXd simi_test = inv(left);
        double erri_test = (matprod(simi_test, left) -
                            Eigen::MatrixXd::Identity(num_vars, num_vars)).cwiseAbs().maxCoeff();
        if (erri_test < erri || (std::isnan(erri) && !std::isnan(erri_test))) {
            simi_new = simi_test;
            erri = erri_test;
        }
    }

    // If SIMI is satisfactory, then update FVAL, CONMAT, CVAL, and the pole position. Otherwise
    // restore SIM and SIMI, and return with INFO = DAMAGING_ROUNDING.
    if (erri <= itol) {
        Eigen::VectorXd fval_new = fval;
        fval_new[jdrop] = f;
        Eigen::MatrixXd conmat_new = conmat;
        conmat_new.col(jdrop) = constr;
        Eigen::VectorXd cval_new = cval;
        cval_new[jdrop] = cstrv;

        // Switch the best vertex to the pole position SIM[:, NUM_VARS] if it is not there already
        std::tie(conmat_new, cval_new, fval_new, sim_new, simi_new, info) =
            updatepole(cpen, conmat_new, cval_new, fval_new, sim_new, simi_new);

        return {sim_new, simi_new, fval_new, conmat_new, cval_new, info};
    } else {
        info = DAMAGING_ROUNDING;
        return {sim, simi, fval, conmat, cval, info};
    }
}

} // namespace prima

#endif // PRIMA_CPP_COBYLA_UPDATE_HPP
