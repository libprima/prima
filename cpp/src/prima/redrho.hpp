#ifndef PRIMA_CPP_REDRHO_HPP
#define PRIMA_CPP_REDRHO_HPP

// This module provides a function that calculates RHO when it needs to be reduced.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <algorithm>
#include <cmath>

namespace prima {

// This function calculates RHO when it needs to be reduced.
// The scheme is shared by UOBYQA, NEWUOA, BOBYQA, LINCOA. For COBYLA, Powell's code reduces RHO by
// 'RHO *= 0.5; if RHO <= 1.5 * RHOEND: RHO = RHOEND' as specified in (11) of the COBYLA
// paper. However, this scheme seems to work better, especially after we introduce DELTA.
inline double redrho(double rho_in, double rhoend) {
    double rho_ratio = rho_in / rhoend;
    double rho;
    if (rho_ratio > 250) rho = 0.1 * rho_in;
    else if (rho_ratio <= 16) rho = rhoend;
    else rho = std::sqrt(rho_ratio) * rhoend;  // rho = sqrt(rho * rhoend)
    return rho;
}

} // namespace prima

#endif // PRIMA_CPP_REDRHO_HPP
