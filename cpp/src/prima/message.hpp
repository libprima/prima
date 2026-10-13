#ifndef PRIMA_CPP_MESSAGE_HPP
#define PRIMA_CPP_MESSAGE_HPP

/*
This module provides some functions that print messages to terminal/files.

Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

N.B.:
1. In case parallelism is desirable (especially during initialization), the functions may
   have to be modified or disabled due to the IO operations.
2. IPRINT indicates the level of verbosity, which increases with the absolute value of IPRINT.
   IPRINT = +/-3 can be expensive due to high IO operations.
*/

#include <Eigen/Core>
#include <cmath>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

#include "infos.hpp"

namespace prima {

// Translate an info code into a human-readable reason string
inline std::string get_info_string(const std::string& solver, int info) {
    std::string reason;
    if (info == FTARGET_ACHIEVED) reason = "the target function value is achieved.";
    else if (info == MAXFUN_REACHED) reason = "the objective function has been evaluated MAXFUN times.";
    else if (info == MAXTR_REACHED) reason = "the maximal number of trust region iterations has been reached.";
    else if (info == SMALL_TR_RADIUS) reason = "the trust region radius reaches its lower bound.";
    else if (info == TRSUBP_FAILED) reason = "a trust region step has failed to reduce the quadratic model.";
    else if (info == NAN_INF_X) reason = "NaN or Inf occurs in x.";
    else if (info == NAN_INF_F) reason = "the objective function returns NaN/+Inf.";
    else if (info == NAN_INF_MODEL) reason = "NaN or Inf occurs in the models.";
    else if (info == DAMAGING_ROUNDING) reason = "rounding errors are becoming damaging.";
    else if (info == NO_SPACE_BETWEEN_BOUNDS) reason = "there is no space between the lower and upper bounds of variable.";
    else if (info == ZERO_LINEAR_CONSTRAINT) reason = "one of the linear constraints has a zero gradient";
    else if (info == CALLBACK_TERMINATE) reason = "the callback function requested termination";
    else reason = "UNKNOWN EXIT FLAG";
    return "Return from " + solver + " because " + reason;
}

// This function prints messages at return
inline void retmsg(const std::string& solver, int info, int iprint, int nf,
                   double f, const Eigen::VectorXd& x, std::optional<double> cstrv = std::nullopt,
                   std::optional<const Eigen::VectorXd*> constr = std::nullopt) {
    if (std::abs(iprint) < 1) return;

    // Decide whether the problem is truly constrained
    bool is_constrained = constr.has_value() || cstrv.has_value();
    // Decide the constraint violation (N.B.: We assume that the constraint is CONSTR >= 0.)
    double cstrv_loc = cstrv.value_or(0);
    if (!cstrv.has_value() && constr.has_value() && constr.value()) {
        double m = 0;
        for (Eigen::Index i = 0; i < constr.value()->size(); ++i)
            m = std::max(m, -(*constr.value())[i]);
        cstrv_loc = m;
    }
    std::string msg = get_info_string(solver, info);

    // Print X in compact form for 2 or fewer elements, otherwise as column
    char buf[1024];
    if (x.size() <= 2) {
        std::snprintf(buf, sizeof(buf), "\nThe corresponding X is: [%g", x[0]);
        for (Eigen::Index i = 1; i < x.size(); ++i)
            std::snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf), ", %g", x[i]);
        std::snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf), "]");
    } else {
        std::snprintf(buf, sizeof(buf), "\nThe corresponding X is:\n");
        for (Eigen::Index i = 0; i < x.size(); ++i)
            std::snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf), "%g ", x[i]);
    }
    msg += buf;

    if (is_constrained)
        std::printf("\nNumber of function values = %d  Least value of F = %g  Constraint violation = %g",
                    nf, f, cstrv_loc);
    else
        std::printf("\nNumber of function values = %d  Least value of F = %g", nf, f);

    std::printf("\n%s\n", msg.c_str());
}

// This subroutine prints messages for each evaluation of the objective function
inline void fmsg(const std::string& solver, const std::string& state, int iprint, int nf,
                 double delta, double f, const Eigen::VectorXd& x,
                 std::optional<double> cstrv = std::nullopt,
                 std::optional<const Eigen::VectorXd*> constr = std::nullopt) {
    if (std::abs(iprint) < 2) return;
    bool is_constrained = constr.has_value() || cstrv.has_value();
    double cstrv_loc = cstrv.value_or(0);
    std::printf("\n%s step with radius = %g", state.c_str(), delta);
    if (is_constrained)
        std::printf("\nNumber of function values = %d  Least value of F = %g  Constraint violation = %g",
                    nf, f, cstrv_loc);
    else
        std::printf("\nNumber of function values = %d  Least value of F = %g", nf, f);
}

// This function prints messages when RHO is updated
inline void rhomsg(const std::string& solver, int iprint, int nf, double f, double rho,
                   const Eigen::VectorXd& x, std::optional<double> cstrv = std::nullopt,
                   std::optional<const Eigen::VectorXd*> constr = std::nullopt,
                   std::optional<double> cpen = std::nullopt) {
    if (std::abs(iprint) < 2) return;
    bool is_constrained = constr.has_value() || cstrv.has_value();
    double cstrv_loc = cstrv.value_or(0);
    if (cpen.has_value())
        std::printf("\nNew RHO = %g  CPEN = %g", rho, cpen.value());
    else
        std::printf("\nNew RHO = %g", rho);
    if (is_constrained)
        std::printf("\nNumber of function values = %d  Least value of F = %g  Constraint violation = %g",
                    nf, f, cstrv_loc);
    else
        std::printf("\nNumber of function values = %d  Least value of F = %g", nf, f);
}

} // namespace prima

#endif // PRIMA_CPP_MESSAGE_HPP
