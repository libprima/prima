#ifndef PRIMA_CPP_EVALUATE_HPP
#define PRIMA_CPP_EVALUATE_HPP

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>

#include "consts.hpp"
#include "linalg.hpp"

// This is a module evaluating the objective/constraint function with Nan/Inf handling.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

namespace prima {

using Calcfc = std::function<std::pair<double, Eigen::VectorXd>(const Eigen::VectorXd&)>;

// This function moderates a decision variable. It replaces NaN by 0 and Inf/-Inf by
// REALMAX/-REALMAX.
inline Eigen::VectorXd moderatex(const Eigen::VectorXd& x) {
    Eigen::VectorXd y = x;
    for (Eigen::Index i = 0; i < y.size(); ++i) {
        if (std::isnan(y[i])) y[i] = 0;
        y[i] = std::max(-REALMAX, std::min(REALMAX, y[i]));
    }
    return y;
}

// This function moderates the function value of a MINIMIZATION problem. It replaces
// NaN and any value above FUNCMAX by FUNCMAX.
// We may moderate huge negative function values as follows, but we decide not to.
//   f = std::max(-FUNCMAX, std::min(FUNCMAX, f));
inline double moderatef(double f) {
    if (std::isnan(f)) f = FUNCMAX;
    return std::max(-REALMAX, std::min(FUNCMAX, f));
}

// This function moderates the constraint value, the constraint demanding this value
// to be NONNEGATIVE. It replaces any value below -CONSTRMAX by -CONSTRMAX, and any
// NaN or value above CONSTRMAX by CONSTRMAX.
inline double moderatec(double c) {
    if (std::isnan(c) || c > CONSTRMAX) c = CONSTRMAX;
    return std::max(-CONSTRMAX, std::min(CONSTRMAX, c));
}

// This function moderates the constraint value, the constraint demanding this value
// to be NONNEGATIVE. It replaces any value below -CONSTRMAX by -CONSTRMAX, and any
// NaN or value above CONSTRMAX by CONSTRMAX.
inline Eigen::VectorXd moderatec(const Eigen::VectorXd& c) {
    Eigen::VectorXd y = c;
    for (Eigen::Index i = 0; i < y.size(); ++i) {
        if (std::isnan(y[i]) || y[i] > CONSTRMAX) y[i] = CONSTRMAX;
        y[i] = std::max(-CONSTRMAX, std::min(CONSTRMAX, y[i]));
    }
    return y;
}

// This function evaluates CALCFC at X, returning the objective function value and the
// constraint value. Nan/Inf are handled by a moderated extreme barrier.
inline std::pair<double, Eigen::VectorXd>
evaluate(const Calcfc& calcfc, const Eigen::VectorXd& x, int m_nlcon,
         const Eigen::MatrixXd* amat, const Eigen::VectorXd* bvec) {

    // Sizes
    int m_lcon = (bvec != nullptr) ? bvec->size() : 0;
    int total_con = m_lcon + m_nlcon;
    Eigen::VectorXd constr(total_con);

    if (amat != nullptr) {
        constr.head(m_lcon) = matprod(x, amat->transpose()) - *bvec;
    }

    double f;
    bool has_nan = false;
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        if (std::isnan(x[i])) { has_nan = true; break; }
    }

    if (has_nan) {
        // Although this should not happen unless there is a bug, we include this case
        // for robustness.
        f = primasum(x);
        constr.tail(m_nlcon).setConstant(f);
    } else {
        auto [fc, cc] = calcfc(moderatex(x));
        f = fc;
        if (m_nlcon > 0) constr.tail(m_nlcon) = cc;

        // Moderated extreme barrier: replace NaN/huge objective or constraint values
        // with a large but finite value. This is naive, and better approaches surely
        // exist.
        f = moderatef(f);
        if (m_nlcon > 0) constr.tail(m_nlcon) = moderatec(constr.tail(m_nlcon));
    }

    return {f, constr};
}

} // namespace prima

#endif // PRIMA_CPP_EVALUATE_HPP
