#ifndef PRIMA_CPP_COBYLA_INITIALIZE_HPP
#define PRIMA_CPP_COBYLA_INITIALIZE_HPP

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

#include "../checkbreak.hpp"
#include "../consts.hpp"
#include "../evaluate.hpp"
#include "../infos.hpp"
#include "../history.hpp"
#include "../linalg.hpp"
#include "../message.hpp"
#include "../selectx.hpp"

namespace prima {

inline std::tuple<Eigen::VectorXi, Eigen::MatrixXd, Eigen::VectorXd,
                  Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, int, int>
initxfc(const Calcfc& calcfc, int iprint, int maxfun,
        const Eigen::VectorXd& constr0, const Eigen::MatrixXd* amat,
        const Eigen::VectorXd* bvec, double ctol, double f0, double ftarget,
        double rhobeg, const Eigen::VectorXd& x0,
        std::vector<Eigen::VectorXd>& xhist, std::vector<double>& fhist,
        std::vector<double>& chist, std::vector<Eigen::VectorXd>& conhist, int maxhist) {
    // This subroutine does the initialization concerning X, function values, and
    // constraints.

    // Local variables
    int num_constraints = constr0.size();
    int m_lcon = bvec ? bvec->size() : 0;
    int m_nlcon = num_constraints - m_lcon;
    int num_vars = x0.size();

    // Initialize info to the default value. At return, a value different from this
    // value will indicate an abnormal return.
    int info = INFO_DEFAULT;

    // Initialize the simplex. It will be revised during the initialization.
    Eigen::MatrixXd sim = Eigen::MatrixXd::Zero(num_vars, num_vars + 1);
    sim.rightCols(1).col(0) = x0;
    for (int i = 0; i < num_vars; ++i) sim(i, i) = rhobeg;

    // Initialize the matrix simi. In most cases simi is overwritten, but not always.
    Eigen::MatrixXd simi = Eigen::MatrixXd::Identity(num_vars, num_vars) / rhobeg;

    // evaluated[j] = true iff the function/constraint of sim.col(j) has been evaluated.
    Eigen::VectorXi evaluated = Eigen::VectorXi::Zero(num_vars + 1);

    // Initialize fval
    Eigen::VectorXd fval = Eigen::VectorXd::Constant(num_vars + 1, REALMAX);
    Eigen::VectorXd cval = Eigen::VectorXd::Constant(num_vars + 1, REALMAX);
    Eigen::MatrixXd conmat = Eigen::MatrixXd::Constant(num_constraints, num_vars + 1, REALMAX);

    for (int k = 0; k < num_vars + 1; ++k) {
        // We will evaluate F corresponding to sim.col(j).
        Eigen::VectorXd x = sim.col(num_vars);
        int j;
        double f;
        Eigen::VectorXd constr;
        if (k == 0) {
            j = num_vars;
            f = f0;
            constr = constr0;
        } else {
            j = k - 1;
            x[j] += rhobeg;
            auto [fv, cv] = evaluate(calcfc, x, m_nlcon, amat, bvec);
            f = fv; constr = cv;
        }
        double cstrv = 0;
        for (int i = 0; i < constr.size(); ++i) cstrv = std::max(cstrv, constr[i]);

        // Print a message about the function/constraint evaluation according to IPRINT.
        fmsg("COBYLA", "Initialization", iprint, k, rhobeg, f, x, cstrv, &constr);

        // Save X, F, CONSTR, CSTRV into the history.
        savehist(maxhist, x, xhist, f, fhist, cstrv, chist, constr, conhist);

        // Save F, CONSTR, and CSTRV to FVAL, CONMAT, and CVAL respectively.
        evaluated[j] = 1;
        fval[j] = f;
        conmat.col(j) = constr;
        cval[j] = cstrv;

        // Check whether to exit.
        int subinfo = checkbreak_con(maxfun, k, cstrv, ctol, f, ftarget, x);
        if (subinfo != INFO_DEFAULT) { info = subinfo; break; }

        // Exchange the new vertex of the initial simplex with the optimal vertex if necessary.
        // This is the ONLY part that is essentially non-parallel.
        if (j < num_vars && fval[j] < fval[num_vars]) {
            std::swap(fval[j], fval[num_vars]);
            std::swap(cval[j], cval[num_vars]);
            conmat.col(j).swap(conmat.col(num_vars));
            sim.col(num_vars) = x;
            // SIM(:, :j+1) is lower triangular
            for (int ii = 0; ii <= j; ++ii) sim(j, ii) = -rhobeg;
        }
    }

    int nf = evaluated.sum();

    if (evaluated.minCoeff() == 1) {
        // Initialize SIMI to the inverse of SIM[:, :num_vars]
        simi = inv(sim.leftCols(num_vars));
    }

    return {evaluated, conmat, cval, sim, simi, fval, nf, info};
}

inline int initfilt(const Eigen::MatrixXd& conmat, double ctol, double cweight,
                    const Eigen::VectorXd& cval, const Eigen::VectorXd& fval,
                    const Eigen::MatrixXd& sim, const Eigen::VectorXi& evaluated,
                    Eigen::VectorXd& cfilt, Eigen::MatrixXd& confilt,
                    Eigen::VectorXd& ffilt, Eigen::MatrixXd& xfilt) {
    // This function initializes the filter (XFILT, etc) that will be used when selecting
    // X at the end of the solver.
    // N.B.:
    // 1. Why not initialize the filters using XHIST, etc? Because the history is empty if
    //    the user chooses not to output it.
    // 2. We decouple INITXFC and INITFILT so that it is easier to parallelize the former
    //    if needed.

    int num_constraints = conmat.rows();
    int num_vars = sim.rows();
    int maxfilt = ffilt.size();
    int nfilt = 0;

    for (int i = 0; i < num_vars + 1; ++i) {
        if (evaluated[i]) {
            Eigen::VectorXd x;
            if (i < num_vars) x = sim.col(i) + sim.col(num_vars);
            else x = sim.col(i);  // i == num_vars, i.e. the last column
            Eigen::VectorXd ci = conmat.col(i);
            nfilt = savefilt(cval[i], ctol, cweight, fval[i], x, nfilt, cfilt, ffilt, xfilt,
                            &ci, &confilt);
        }
    }

    return nfilt;
}

} // namespace prima

#endif // PRIMA_CPP_COBYLA_INITIALIZE_HPP
