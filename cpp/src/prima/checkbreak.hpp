#ifndef PRIMA_CPP_CHECKBREAK_HPP
#define PRIMA_CPP_CHECKBREAK_HPP

// This module checks whether to break out of the solver loop.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <Eigen/Core>
#include <algorithm>
#include <cmath>

#include "infos.hpp"
#include "consts.hpp"

namespace prima {

// This function checks whether to break out of the solver loop in the unconstrained case.
inline int checkbreak_unc(int maxfun, int nf, double f, double ftarget, const Eigen::VectorXd& x) {
    int info = INFO_DEFAULT;

    // Although X should not contain NaN unless there is a bug, we include the following for security.
    // X can be Inf, as finite + finite can be Inf numerically.
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        if (std::isnan(x[i]) || std::isinf(x[i])) { info = NAN_INF_X; break; }
    }

    // Although NAN_INF_F should not happen unless there is a bug, we include the following for security.
    if (std::isnan(f) || std::isinf(f)) info = NAN_INF_F;

    if (f <= ftarget) info = FTARGET_ACHIEVED;
    if (nf >= maxfun) info = MAXFUN_REACHED;

    return info;
}

// This function checks whether to break out of the solver loop in the constrained case.
inline int checkbreak_con(int maxfun, int nf, double cstrv, double ctol,
                          double f, double ftarget, const Eigen::VectorXd& x) {
    int info = INFO_DEFAULT;

    // Although X should not contain NaN unless there is a bug, we include the following for security.
    // X can be Inf, as finite + finite can be Inf numerically.
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        if (std::isnan(x[i]) || std::isinf(x[i])) { info = NAN_INF_X; break; }
    }

    // Although NAN_INF_F should not happen unless there is a bug, we include the following for security.
    if (std::isnan(f) || std::isinf(f) || std::isnan(cstrv) || std::isinf(cstrv))
        info = NAN_INF_F;

    if (cstrv <= ctol && f <= ftarget) info = FTARGET_ACHIEVED;
    if (nf >= maxfun) info = MAXFUN_REACHED;

    return info;
}

} // namespace prima

#endif // PRIMA_CPP_CHECKBREAK_HPP
