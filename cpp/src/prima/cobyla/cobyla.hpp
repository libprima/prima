#ifndef PRIMA_CPP_COBYLA_COBYLA_HPP
#define PRIMA_CPP_COBYLA_COBYLA_HPP

// This module provides Powell's COBYLA algorithm.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
//
// N.B.:
//
// 1. The modern-Fortran reference implementation in PRIMA contains bug fixes and improvements over the
// original Fortran 77 implementation by Powell. Consequently, the PRIMA implementation behaves differently
// from the original Fortran 77 implementation by Powell. Therefore, it is important to point out that
// you are using PRIMA rather than the original solvers if you want your results to be reproducible.
//
// 2. Compared to Powell's Fortran 77 implementation, the modern-Fortran implementation and hence any
// faithful translation like this one generally produce better solutions with fewer function evaluations,
// making them preferable for applications with expensive function evaluations. However, if function
// evaluations are not the dominant cost in your application, the Fortran 77 solvers are likely to be
// faster, as they are more efficient in terms of memory usage and flops thanks to the careful and
// ingenious (but unmaintained and unmaintainable) implementation by Powell.
//
// See the PRIMA documentation (www.libprima.net) for more information.

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>
#include <tuple>
#include <vector>

#include "../consts.hpp"
#include "../evaluate.hpp"
#include "../linalg.hpp"
#include "../preproc.hpp"
#include "../present.hpp"

#include "cobylb.hpp"

namespace prima {

struct COBYLAResult {
    Eigen::VectorXd x;
    double f;
    Eigen::VectorXd constr;
    double cstrv;
    int nf;
    std::vector<Eigen::VectorXd> xhist;
    std::vector<double> fhist;
    std::vector<double> chist;
    std::vector<Eigen::VectorXd> conhist;
    int info;
};

// This subroutine wraps the linear and bound constraints into a single constraint:
//     AMAT*X <= BVEC.
//
// N.B.:
//
// LINCOA normalizes the linear constraints so that each constraint has a gradient
// of norm 1. However, COBYLA does not do this.
inline std::tuple<std::optional<Eigen::MatrixXd>, std::optional<Eigen::VectorXd>>
get_lincon(const std::optional<Eigen::MatrixXd>& Aeq,
           const std::optional<Eigen::MatrixXd>& Aineq,
           const std::optional<Eigen::VectorXd>& beq,
           const std::optional<Eigen::VectorXd>& bineq,
           const std::optional<Eigen::VectorXd>& xl,
           const std::optional<Eigen::VectorXd>& xu) {
    // Sizes
    int num_vars = 0;
    if (Aeq.has_value()) num_vars = Aeq->cols();
    else if (Aineq.has_value()) num_vars = Aineq->cols();
    else if (xl.has_value()) num_vars = xl->size();
    else if (xu.has_value()) num_vars = xu->size();
    else return {std::nullopt, std::nullopt};

    //====================//
    // Calculation starts //
    //====================//

    // Define the indices of the nontrivial bound constraints.
    std::vector<int> ixl;
    std::vector<int> ixu;

    if (xl.has_value()) {
        for (int i = 0; i < xl->size(); ++i)
            if ((*xl)[i] > -BOUNDMAX) ixl.push_back(i);
    }
    if (xu.has_value()) {
        for (int i = 0; i < xu->size(); ++i)
            if ((*xu)[i] < BOUNDMAX) ixu.push_back(i);
    }

    // Wrap the linear constraints.
    // The bound constraint XL <= X <= XU is handled as two constraints:
    // -X <= -XL, X <= XU.
    // The equality constraint Aeq*X = Beq is handled as two constraints:
    // -Aeq*X <= -Beq, Aeq*X <= Beq.
    // N.B.:
    // 1. The treatment of the equality constraints is naive. One may choose to
    //    eliminate them instead.
    Eigen::MatrixXd idmat = Eigen::MatrixXd::Identity(num_vars, num_vars);

    int mxl = ixl.size();
    int mxu = ixu.size();
    int meq = Aeq.has_value() ? Aeq->rows() : 0;
    int mineq = Aineq.has_value() ? Aineq->rows() : 0;
    int total_rows = mxl + mxu + 2 * meq + mineq;

    if (total_rows == 0) return {std::nullopt, std::nullopt};

    Eigen::MatrixXd amat(total_rows, num_vars);
    Eigen::VectorXd bvec(total_rows);

    int row = 0;
    for (int i : ixl) {
        amat.row(row) = -idmat.row(i);
        bvec[row] = -(*xl)[i];
        row++;
    }
    for (int i : ixu) {
        amat.row(row) = idmat.row(i);
        bvec[row] = (*xu)[i];
        row++;
    }
    if (Aeq.has_value()) {
        for (int i = 0; i < Aeq->rows(); ++i) {
            amat.row(row) = -Aeq->row(i);
            bvec[row] = -(*beq)[i];
            row++;
        }
        for (int i = 0; i < Aeq->rows(); ++i) {
            amat.row(row) = Aeq->row(i);
            bvec[row] = (*beq)[i];
            row++;
        }
    }
    if (Aineq.has_value()) {
        for (int i = 0; i < Aineq->rows(); ++i) {
            amat.row(row) = Aineq->row(i);
            bvec[row] = (*bineq)[i];
            row++;
        }
    }

    //==================//
    // Calculation ends //
    //==================//

    return {amat, bvec};
}

// Among all the arguments, only CALCFC, M_NLCON, and X are obligatory. The others are
// OPTIONAL and you can neglect them unless you are familiar with the algorithm. Any
// unspecified optional input will take the default value detailed below. For
// instance, we may invoke the solver as follows.
//
//   // First define CALCFC, M_NLCON, and X, and then do the following.
//   result = cobyla(calcfc, m_nlcon, x);
//
// or
//
//   // First define CALCFC, M_NLCON, X, Aineq, and Bineq, and then do the following.
//   result = cobyla(calcfc, m_nlcon, x, Aineq, bineq, rhobeg=1.0e0,
//       rhoend=1.0e-6);
//
// ####################################################################################
// # IMPORTANT NOTICE: The user must set M_NLCON correctly to the number of nonlinear
// # constraints, namely the size of NLCONSTR introduced below. Set it to 0 if there
// # is no nonlinear constraint.
// ####################################################################################
//
// See examples/cobyla/cobyla_example.py for a concrete example.
//
// A detailed introduction to the arguments is as follows.
//
// ####################################################################################
// # INPUTS
// ####################################################################################
//
// CALCFC
//   Input, function.
//   f, nlconstr = CALCFC(X) should evaluate the objective function and nonlinear
//   constraints at the given vector X; it should return a tuple consisting of the
//   objective function value and the nonlinear constraint value. It must be provided
//   by the user, and its definition must conform to the following interface:
//   //-------------------------------------------------------------------------//
//    auto [f, nlconstr] = calcfc(x);
//   //-------------------------------------------------------------------------//
//
// M_NLCON
//   Input, scalar.
//   M_NLCON must be set to the number of nonlinear constraints, namely the size of
//   NLCONSTR(X).
//   N.B.:
//   1. Why don't we define M_NLCON as optional and default it to 0 when it is absent?
//   This is because we need to allocate memory for CONSTR_LOC using M_NLCON. To
//   ensure that the size of CONSTR_LOC is correct, we require the user to specify
//   M_NLCON explicitly.
//
// X
//   Input, vector.
//   As an input, X should be an N-dimensional vector that contains the starting
//   point, N being the dimension of the problem.
//
// Aineq, Bineq
//   Input, matrix of size [Mineq, N] and vector of size Mineq unless they are both
//   empty, default: None and None.
//   Aineq and Bineq represent the linear inequality constraints: Aineq*X <= Bineq.
//
// Aeq, Beq
//   Input, matrix of size [Meq, N] and vector of size Meq unless they are both
//   empty, default: None and None.
//   Aeq and Beq represent the linear equality constraints: Aeq*X = Beq.
//
// XL, XU
//   Input, vectors of size N unless they are both None, default: None and None.
//   XL is the lower bound for X. If XL is None, X has no
//   lower bound. Any entry of XL that is NaN or below -BOUNDMAX will be taken as
//   -BOUNDMAX, which effectively means there is no lower bound for the corresponding
//   entry of X. The value of BOUNDMAX is 0.25*HUGE(X), which is about 8.6E37 for
//   single precision and 4.5E307 for double precision. XU is similar.
//
// F0
//   Input, scalar.
//   F0, if present, should be set to the objective function value of the starting X.
//
// NLCONSTR0
//   Input, vector.
//   NLCONSTR0, if present, should be set to the nonlinear constraint value at the
//   starting X; in addition, SIZE(NLCONSTR0) must be M_NLCON, or the solver will
//   abort.
//
// RHOBEG, RHOEND
//   Inputs, scalars, default: RHOBEG = 1, RHOEND = 10^-6. RHOBEG and RHOEND must be
//   set to the initial and final values of a trust-region radius, both being positive
//   and RHOEND <= RHOBEG. Typically RHOBEG should be about one tenth of the greatest
//   expected change to a variable, and RHOEND should indicate the accuracy that is
//   required in the final values of the variables.
//
// FTARGET
//   Input, scalar, default: -Inf.
//   FTARGET is the target function value. The algorithm will terminate when a
//   feasible point with a function value <= FTARGET is found.
//
// CTOL
//   Input, scalar, default: sqrt(machine epsilon).
//   CTOL is the tolerance of constraint violation. X is considered feasible if
//   CSTRV(X) <= CTOL.
//   N.B.:
//     1. CTOL is absolute, not relative.
//     2. CTOL is used only when selecting the returned X. It does not affect the
//        iterations of the algorithm.
//
// CWEIGHT
//   Input, scalar, default: CWEIGHT_DFT defined in common/consts.hpp.
//   CWEIGHT is the weight that the constraint violation takes in the selection of the
//   returned X.
//
// MAXFUN
//   Input, integer scalar, default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in
//   common/consts.hpp. MAXFUN is the maximal number of calls of CALCFC.
//
// IPRINT
//   Input, integer scalar, default: 0.
//   The value of IPRINT should be set to 0, 1, -1, 2, -2, 3, or -3, which controls
//   how much information will be printed during the computation:
//   0: there will be no printing;
//   1: a message will be printed to the screen at the return, showing the best vector
//      of variables found and its objective function value;
//   2: in addition to 1, each new value of RHO is printed to the screen, with the
//      best vector of variables so far and its objective function value; each new
//      value of CPEN is also printed;
//   3: in addition to 2, each function evaluation with its variables will be printed
//      to the screen; -1, -2, -3: the same information as 1, 2, 3 will be printed,
//      not to the screen but to a file named COBYLA_output.txt; the file will be
//      created if it does not exist; the new output will be appended to the end of
//      this file if it already exists.
//   Note that IPRINT = +/-3 can be costly in terms of time and/or space.
//
// ETA1, ETA2, GAMMA1, GAMMA2
//   Input, scalars, default: ETA1 = 0.1, ETA2 = 0.7, GAMMA1 = 0.5, and GAMMA2 = 2.
//   ETA1, ETA2, GAMMA1, and GAMMA2 are parameters in the updating scheme of the
//   trust-region radius detailed in the subroutine TRRAD in trustregion.hpp. Roughly
//   speaking, the trust-region radius is contracted by a factor of GAMMA1 when the
//   reduction ratio is below ETA1, and enlarged by a factor of GAMMA2 when the
//   reduction ratio is above ETA2. It is required that 0 < ETA1 <= ETA2 < 1 and
//   0 < GAMMA1 < 1 < GAMMA2. Normally, ETA1 <= 0.25. It is NOT advised to set
//   ETA1 >= 0.5.
//
// MAXFILT
//   Input, scalar.
//   MAXFILT is a nonnegative integer indicating the maximal length of the filter used
//   for selecting the returned solution; default: MAXFILT_DFT (a value lower than
//   MIN_MAXFILT is not recommended);
//   see common/consts.hpp for the definitions of MAXFILT_DFT and MIN_MAXFILT.
//
// CALLBACK
//   Input, function to report progress and optionally request termination.
//
//
// ####################################################################################
// # OUTPUTS
// ####################################################################################
//
// The output is a single data structure, COBYLAResult, with the following fields:
//
// X
//   Output, vector.
//   As an output, X will be set to an approximate minimizer.
//
// F
//   Output, scalar.
//   F will be set to the objective function value of X at exit.
//
// CONSTR
//   Output, vector.
//   CONSTR will be set to the constraint value of X at exit.
//
// CSTRV
//   Output, scalar.
//   CSTRV will be set to the constraint violation of X at exit, i.e.,
//   max([0, XL - X, X - XU, Aineq*X - Bineq, ABS(Aeq*X - Beq), NLCONSTR(X)]).
//
// NF
//   Output, scalar.
//   NF will be set to the number of calls of CALCFC at exit.
//
// XHIST, FHIST, CHIST, CONHIST, MAXHIST
//   XHIST: Output, rank 2 array;
//   FHIST: Output, rank 1 array;
//   CHIST: Output, rank 1 array;
//   CONHIST: Output, rank 2 array;
//   MAXHIST: Input, scalar, default: MAXFUN
//   XHIST, if present, will output the history of iterates; FHIST, if present, will
//   output the history function values; CHIST, if present, will output the history of
//   constraint violations; CONHIST, if present, will output the history of constraint
//   values; MAXHIST should be a nonnegative integer, and XHIST/FHIST/CHIST/CONHIST
//   will output only the history of the last MAXHIST iterations.
//   Therefore, MAXHIST= 0 means XHIST/FHIST/CONHIST/CHIST will output
//   nothing, while setting MAXHIST = MAXFUN requests XHIST/FHIST/CHIST/CONHIST to
//   output all the history.
//
//   IMPORTANT NOTICE:
//   Setting MAXHIST to a large value can be costly in terms of memory for large
//   problems.
//   MAXHIST will be reset to a smaller value if the memory needed exceeds MAXHISTMEM
//   defined in common/consts.hpp.
//   Use *HIST with caution!!! (N.B.: the algorithm is NOT designed for large
//   problems).
//
// INFO
//   Output, scalar.
//   INFO is the exit flag. It will be set to one of the following values defined in
//   common/infos.hpp:
//   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
//   FTARGET_ACHIEVED: the target function value is reached;
//   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
//   MAXTR_REACHED: the trust region iteration has been performed MAXTR times (MAXTR = 2*MAXFUN);
//   NAN_INF_X: NaN or Inf occurs in X;
//   DAMAGING_ROUNDING: rounding errors are becoming damaging.
//   //--------------------------------------------------------------------------//
//   The following case(s) should NEVER occur unless there is a bug.
//   NAN_INF_F: the objective function returns NaN or +Inf;
//   NAN_INF_MODEL: NaN or Inf occurs in the model;
//   TRSUBP_FAILED: a trust region step failed to reduce the model
//   //--------------------------------------------------------------------------//
inline COBYLAResult
cobyla(const Calcfc& calcfc, int m_nlcon, const Eigen::VectorXd& x,
       std::optional<Eigen::MatrixXd> Aineq = std::nullopt,
       std::optional<Eigen::VectorXd> bineq = std::nullopt,
       std::optional<Eigen::MatrixXd> Aeq = std::nullopt,
       std::optional<Eigen::VectorXd> beq = std::nullopt,
       std::optional<Eigen::VectorXd> xl = std::nullopt,
       std::optional<Eigen::VectorXd> xu = std::nullopt,
       std::optional<double> f0 = std::nullopt,
       std::optional<Eigen::VectorXd> nlconstr0 = std::nullopt,
       std::optional<double> rhobeg = std::nullopt,
       std::optional<double> rhoend = std::nullopt,
       double ftarget = FTARGET_DEFAULT,
       double ctol = CTOL_DEFAULT,
       double cweight = CWEIGHT_DEFAULT,
       std::optional<int> maxfun = std::nullopt,
       int iprint = IPRINT_DEFAULT,
       std::optional<double> eta1 = std::nullopt,
       std::optional<double> eta2 = std::nullopt,
       double gamma1 = GAMMA1_DEFAULT,
       double gamma2 = GAMMA2_DEFAULT,
       std::optional<int> maxhist = std::nullopt,
       int maxfilt = 2000,
       std::function<bool(const Eigen::VectorXd&, double, int, int, double, const Eigen::VectorXd&)> callback = nullptr) {

    // Local variables
    std::string solver = "COBYLA";
    int num_vars = x.size();

    // Sizes
    int mineq = Aineq.has_value() ? Aineq->rows() : 0;
    int meq = Aeq.has_value() ? Aeq->rows() : 0;
    int mxl = 0, mxu = 0;
    if (xl.has_value())
        for (int i = 0; i < xl->size(); ++i)
            if ((*xl)[i] > -BOUNDMAX) mxl++;
    if (xu.has_value())
        for (int i = 0; i < xu->size(); ++i)
            if ((*xu)[i] < BOUNDMAX) mxu++;
    int mmm = mxl + mxu + 2 * meq + mineq + m_nlcon;

    // Wrap the linear and bound constraints into a single constraint: AMAT@X <= BVEC.
    auto [amat_opt, bvec_opt] = get_lincon(Aeq, Aineq, beq, bineq, xl, xu);

    // Create constraint vector
    Eigen::VectorXd constr = Eigen::VectorXd::Zero(mmm);

    // Set [F_LOC, CONSTR_LOC] to [F(X0), CONSTR(X0)] after evaluating the latter if
    // needed. In this way, COBYLB only needs one interface.
    // N.B.: Due to the preconditions above, there are two possibilities for F0 and
    // NLCONSTR0.
    // If NLCONSTR0 is present, then F0 must be present, and we assume that F(X0) = F0
    // even if F0 is NaN.
    // If NLCONSTR0 is absent, then F0 must be either absent or NaN, both of which will
    // be interpreted as F(X0) is not provided and we have to evaluate F(X0) and
    // NLCONSTR(X0) now.
    Eigen::VectorXd constr_loc;
    double f_loc;
    if (f0.has_value() && nlconstr0.has_value()) {
        bool x_finite = true;
        for (int i = 0; i < x.size(); ++i)
            if (!std::isfinite(x[i])) { x_finite = false; break; }
        if (x_finite) {
            f_loc = moderatef(f0.value());
            if (amat_opt.has_value()) {
                constr.head(mmm - m_nlcon) = moderatec(matprod(amat_opt.value(), x) - bvec_opt.value());
            }
            constr.tail(m_nlcon) = moderatec(nlconstr0.value());
        } else {
            Eigen::VectorXd xm = moderatex(x);
            std::tie(f_loc, constr_loc) = evaluate(calcfc, xm, m_nlcon,
                                                      amat_opt.has_value() ? &amat_opt.value() : nullptr,
                                                      bvec_opt.has_value() ? &bvec_opt.value() : nullptr);
            constr = constr_loc;
        }
    } else {
        Eigen::VectorXd xm = moderatex(x);
        std::tie(f_loc, constr_loc) = evaluate(calcfc, xm, m_nlcon,
                                                  amat_opt.has_value() ? &amat_opt.value() : nullptr,
                                                  bvec_opt.has_value() ? &bvec_opt.value() : nullptr);
        constr = constr_loc;
    }
    // N.B.: Do NOT call FMSG, SAVEHIST, or SAVEFILT for the function/constraint evaluation at X0.
    // They will be called during the initialization, which will read the function/constraint at X0.
    double cstrv = 0;
    for (int i = 0; i < constr.size(); ++i)
        cstrv = std::max(cstrv, constr[i]);

    // If RHOBEG is present, use it; otherwise, RHOBEG takes the default value for
    // RHOBEG, taking the value of RHOEND into account. Note that RHOEND is considered
    // only if it is present and it is VALID (i.e., finite and positive). The other
    // inputs are read similarly.
    double rhobeg_loc;
    if (rhobeg.has_value()) {
        rhobeg_loc = rhobeg.value();
    } else if (rhoend.has_value() && std::isfinite(rhoend.value()) && rhoend.value() > 0) {
        rhobeg_loc = std::max(10.0 * rhoend.value(), RHOBEG_DEFAULT);
    } else {
        rhobeg_loc = RHOBEG_DEFAULT;
    }

    double rhoend_loc;
    if (rhoend.has_value()) {
        rhoend_loc = rhoend.value();
    } else if (rhobeg_loc > 0) {
        rhoend_loc = std::max(EPS, std::min(RHOEND_DEFAULT / RHOBEG_DEFAULT * rhobeg_loc, RHOEND_DEFAULT));
    } else {
        rhoend_loc = RHOEND_DEFAULT;
    }

    int maxfun_loc = maxfun.has_value() ? maxfun.value() : MAXFUN_DIM_DEFAULT * num_vars;

    double eta1_loc;
    if (eta1.has_value()) {
        eta1_loc = eta1.value();
    } else if (eta2.has_value() && 0 < eta2.value() && eta2.value() < 1) {
        eta1_loc = std::max(EPS, eta2.value() / 7.0);
    } else {
        eta1_loc = ETA1_DEFAULT;
    }

    double eta2_loc;
    if (eta2.has_value()) {
        eta2_loc = eta2.value();
    } else if (0 < eta1_loc && eta1_loc < 1) {
        eta2_loc = (eta1_loc + 2.0) / 3.0;
    } else {
        eta2_loc = ETA2_DEFAULT;
    }

    int maxhist_loc = maxhist.has_value() ? maxhist.value()
        : std::max({maxfun_loc, num_vars + 2, MAXFUN_DIM_DEFAULT * num_vars});

    // Preprocess the inputs in case some of them are invalid. It does nothing if all
    // inputs are valid.
    PreprocResult pp = preproc(solver, num_vars, iprint, maxfun_loc, maxhist_loc,
                                ftarget, rhobeg_loc, rhoend_loc,
                                mmm, std::nullopt, maxfilt, ctol, cweight,
                                eta1_loc, eta2_loc, gamma1, gamma2, (mmm > 0));
    // npt is unused in COBYLA
    // _x0 is unused in COBYLA

    // Further revise MAXHIST according to MAXHISTMEM, and allocate memory for the history.
    // In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
    // CHIST = NaN(1, MAXFUN), CONHIST = NaN(M, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
    // if they are requested; replace MAXFUN with 0 for the history that is not requested.

    // call cobylb, which performs the real calculations
    auto [x_out, f_out, constr_out, cstrv_out, nf_out,
          xhist_out, fhist_out, chist_out, conhist_out, info_out] =
        cobylb(calcfc, pp.iprint, pp.maxfilt, pp.maxfun,
               amat_opt.has_value() ? &amat_opt.value() : nullptr,
               bvec_opt.has_value() ? &bvec_opt.value() : nullptr,
               pp.ctol, pp.cweight, pp.eta1, pp.eta2,
               pp.ftarget, pp.gamma1, pp.gamma2,
               pp.rhobeg, pp.rhoend,
               constr, f_loc, x, pp.maxhist, callback);

    Eigen::VectorXd final_constr = constr_out.tail(m_nlcon);

    return {x_out, f_out, final_constr, cstrv_out, nf_out,
            xhist_out, fhist_out, chist_out, conhist_out, info_out};
}

} // namespace prima

#endif // PRIMA_CPP_COBYLA_COBYLA_HPP
