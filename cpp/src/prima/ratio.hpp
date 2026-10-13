#ifndef PRIMA_CPP_RATIO_HPP
#define PRIMA_CPP_RATIO_HPP

// This module calculates the reduction ratio for trust-region methods.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <algorithm>
#include <cmath>
#include <limits>

#include "consts.hpp"

namespace prima {

// This function evaluates the reduction ratio of a trust-region step, handling inf/nan properly.
inline double redrat(double ared, double pred, double rshrink) {
    if (std::isnan(ared)) {
        // This should not happen in unconstrained problems due to the moderated extreme barrier.
        return -REALMAX;
    }
    if (std::isnan(pred) || pred <= 0) {
        // The trust-region subproblem solver fails in this rare case. Instead of terminating as
        // Powell's original code does, we set ratio as follows so that the solver may continue
        // to progress.
        if (ared > 0) {
            // The trial point will be accepted, but the trust-region radius will be shrunk if
            // rshrink > 0.
            return rshrink / 2;
        }
        // Set the ratio to a large negative number to signify a bad trust-region step, so that the
        // solver will check whether to take a geometry step or reduce rho.
        return -REALMAX;
    }
    if (std::isinf(pred) && std::isinf(ared)) return 1;       // ared/pred = NaN if calculated directly
    if (std::isinf(pred) && std::isinf(-ared)) return -REALMAX;  // ared/pred = NaN if calculated directly
    return ared / pred;
}

} // namespace prima

#endif // PRIMA_CPP_RATIO_HPP
