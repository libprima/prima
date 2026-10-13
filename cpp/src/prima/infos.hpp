#ifndef PRIMA_CPP_INFOS_HPP
#define PRIMA_CPP_INFOS_HPP

// This is a module defining exit flags.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

namespace prima {

constexpr int INFO_DEFAULT = 0;
constexpr int SMALL_TR_RADIUS = 0;
constexpr int FTARGET_ACHIEVED = 1;
constexpr int TRSUBP_FAILED = 2;
constexpr int MAXFUN_REACHED = 3;
constexpr int FIXED_SUCCESS = 13;
constexpr int MAXTR_REACHED = 20;
constexpr int NAN_INF_X = -1;
constexpr int NAN_INF_F = -2;
constexpr int NAN_INF_MODEL = -3;
constexpr int NO_SPACE_BETWEEN_BOUNDS = 6;
constexpr int DAMAGING_ROUNDING = 7;
constexpr int ZERO_LINEAR_CONSTRAINT = 8;
constexpr int CALLBACK_TERMINATE = 30;

// Stop-codes.
// The following codes are used by ERROR STOP as stop-codes, which should be default integers.
constexpr int INVALID_INPUT = 100;
constexpr int ASSERTION_FAILS = 101;
constexpr int VALIDATION_FAILS = 102;
constexpr int MEMORY_ALLOCATION_FAILS = 103;

} // namespace prima

#endif // PRIMA_CPP_INFOS_HPP
