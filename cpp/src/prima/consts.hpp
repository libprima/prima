#ifndef PRIMA_CPP_CONSTS_HPP
#define PRIMA_CPP_CONSTS_HPP

/*
This is a module defining some constants.

Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
*/

#include <cmath>
#include <limits>

namespace prima {

inline constexpr bool DEBUGGING = false;

// Smallest positive normalized double (denorm_min gives smallest subnormal)
inline constexpr double REALMIN = std::numeric_limits<double>::denorm_min();
inline constexpr double REALMAX = std::numeric_limits<double>::max();
inline constexpr double FUNCMAX = 1.0e30;
inline constexpr double CONSTRMAX = FUNCMAX;
inline double EPS = std::numeric_limits<double>::epsilon();

// Any bound with an absolute value at least BOUNDMAX is considered as no bound.
inline constexpr double BOUNDMAX = REALMAX / 4.0;

// Some default values
inline constexpr double RHOBEG_DEFAULT = 1.0;
inline constexpr double RHOEND_DEFAULT = 1e-6;
inline constexpr double FTARGET_DEFAULT = -REALMAX;
inline const double CTOL_DEFAULT = std::sqrt(std::numeric_limits<double>::epsilon());
inline constexpr double CWEIGHT_DEFAULT = 1e8;
inline constexpr double ETA1_DEFAULT = 0.1;
inline constexpr double ETA2_DEFAULT = 0.7;
inline constexpr double GAMMA1_DEFAULT = 0.5;
inline constexpr double GAMMA2_DEFAULT = 2.0;
inline constexpr int IPRINT_DEFAULT = 0;
inline constexpr int MAXFUN_DIM_DEFAULT = 500;

// 1MB > 10^5 * REAL64. 100 can be too small.
inline constexpr int PRIMA_MAX_HIST_MEM_MB = 300;

// Maximal amount of memory (Byte) allowed for XHIST, FHIST, CONHIST, CHIST, and the filters.
inline constexpr double MHM = PRIMA_MAX_HIST_MEM_MB * 1.0e6;

// Make sure that MAXHISTMEM does not exceed HUGE(0) to avoid overflow and memory errors.
inline constexpr int MAXHISTMEM = static_cast<int>(std::min(MHM, static_cast<double>(std::numeric_limits<int>::max())));

// Maximal length of the filter used in constrained solvers. Should be positive; < 200 is not recommended.
inline constexpr int MIN_MAXFILT = 200;
inline constexpr int MAXFILT_DEFAULT = 10 * MIN_MAXFILT;

} // namespace prima

#endif // PRIMA_CPP_CONSTS_HPP
