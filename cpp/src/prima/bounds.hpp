#ifndef PRIMA_CPP_BOUNDS_HPP
#define PRIMA_CPP_BOUNDS_HPP

/*
This module handles bound constraints for optimization variables.

bounds can either be an object with the properties lb and ub, or a list of tuples
indicating a lower bound and an upper bound for each variable. If the list contains
fewer entries than the length of x0, the remaining entries will be generated as -/+ infinity.

Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
*/

#include <Eigen/Core>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace prima {

// Container for lower and upper bounds
struct Bounds {
    Eigen::VectorXd lb;
    Eigen::VectorXd ub;

    Bounds() = default;

    Bounds(const Eigen::VectorXd& lb_, const Eigen::VectorXd& ub_)
        : lb(lb_), ub(ub_) {}
};

/*
Process bounds for optimization:
- If bounds is nullptr, returns [-inf, inf] for all variables.
- Truncates bounds if longer than lenx0.
- Pads bounds with -/+inf if shorter than lenx0.
- Checks for infeasible bounds (lb > ub).

Returns a pair of (lb, ub) vectors of length lenx0.
*/
inline std::pair<Eigen::VectorXd, Eigen::VectorXd>
process_bounds(const Bounds* bounds, int lenx0) {
    using namespace Eigen;

    if (bounds == nullptr) {
        VectorXd lb = VectorXd::Constant(lenx0, -std::numeric_limits<double>::infinity());
        VectorXd ub = VectorXd::Constant(lenx0, std::numeric_limits<double>::infinity());
        return {lb, ub};
    }

    VectorXd lb = bounds->lb;
    VectorXd ub = bounds->ub;

    // Truncate if bounds are longer than lenx0
    if (lb.size() > lenx0) lb = lb.head(lenx0);
    if (ub.size() > lenx0) ub = ub.head(lenx0);

    // Pad if bounds are shorter than lenx0
    if (lb.size() < lenx0) {
        VectorXd pad = VectorXd::Constant(lenx0 - lb.size(), -std::numeric_limits<double>::infinity());
        VectorXd tmp(lb.size() + pad.size());
        tmp << lb, pad;
        lb = tmp;
    }
    if (ub.size() < lenx0) {
        VectorXd pad = VectorXd::Constant(lenx0 - ub.size(), std::numeric_limits<double>::infinity());
        VectorXd tmp(ub.size() + pad.size());
        tmp << ub, pad;
        ub = tmp;
    }

    // Check the infeasibility of the bounds
    for (int i = 0; i < lenx0; ++i) {
        if (lb(i) > ub(i)) {
            throw std::invalid_argument(
                "Some of the provided bounds are infeasible."
            );
        }
    }

    return {lb, ub};
}

} // namespace prima

#endif // PRIMA_CPP_BOUNDS_HPP
