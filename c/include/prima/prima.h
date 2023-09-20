// Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).

#ifndef PRIMA_H
#define PRIMA_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
# define PRIMAC_HELPER_DLL_IMPORT __declspec(dllimport)
# define PRIMAC_HELPER_DLL_EXPORT __declspec(dllexport)
#else
# define PRIMAC_HELPER_DLL_IMPORT
# define PRIMAC_HELPER_DLL_EXPORT
#endif

#ifdef PRIMAC_STATIC
# define PRIMAC_API
#elif defined primac_EXPORTS
# define PRIMAC_API PRIMAC_HELPER_DLL_EXPORT
#else
# define PRIMAC_API PRIMAC_HELPER_DLL_IMPORT
#endif

/*
 * Verbosity level
 */
typedef enum {
  PRIMA_MSG_NONE = 0, /* No messages */
  PRIMA_MSG_EXIT = 1, /* Exit reasons */
  PRIMA_MSG_RHO = 2, /* Rho changes */
  PRIMA_MSG_FEVL = 3, /* The object/constraint functions get evaluated */
} prima_message;

/*
 * Possible return values
 */
typedef enum
{
  PRIMA_SMALL_TR_RADIUS = 0,
  PRIMA_FTARGET_ACHIEVED = 1,
  PRIMA_TRSUBP_FAILED = 2,
  PRIMA_MAXFUN_REACHED = 3,
  PRIMA_MAXTR_REACHED = 20,
  PRIMA_NAN_INF_X = -1,
  PRIMA_NAN_INF_F = -2,
  PRIMA_NAN_INF_MODEL = -3,
  PRIMA_NO_SPACE_BETWEEN_BOUNDS = 6,
  PRIMA_DAMAGING_ROUNDING = 7,
  PRIMA_ZERO_LINEAR_CONSTRAINT = 8,
  PRIMA_INVALID_INPUT = 100,
  PRIMA_ASSERTION_FAILS = 101,
  PRIMA_VALIDATION_FAILS = 102,
  PRIMA_MEMORY_ALLOCATION_FAILS = 103,
} prima_rc;

/*
 * Return code string
 */
PRIMAC_API
const char *prima_get_rc_string(int rc);

/*
 * A function as required by solvers
 *
 * x     : on input, then vector of variables (should not be modified)
 * f     : on output, the value of the function
 *         a NaN value can be passed to signal an evaluation error
 * constr : on output, the value of the constraints (of size m_nlcon)
 *          NaN values can be passed to signal an evaluation error
 *          only for cobyla
*/
typedef void (*prima_obj)(const double x[], double *f);
typedef void (*prima_objcon)(const double x[], double *f, double constr[]);

/*
 * calfun    : function to minimize (see prima_obj)
 * n         : number of variables (>=0)
 * x         : on input, initial estimate
 *             on output, the solution
 * f         : objective value (output)
 * nf        : number of objective function calls (output)
 * rhobeg    : a reasonable initial change to the variables
 * rhoend    : required accuracy for the variables
 * ftarget   : target function value; optimization stops when f <= ftarget for a feasible point,
 *             can be set to -INFINITY to disable
 * maxfun    : maximum number of function evaluations
 * npt       : number of points in the interpolation set, n+2<=npt<=(n+1)(n+2)/2, recommended: 2*n+1
 * iprint    : verbosity level, see the prima_message enum
 * m_nlcon   : number of non-linear constraints (>=0)
 * calcfc    : function to minimize and constraints (see prima_objcon)
 * cstrv     : constraint violation (output)
 * nlconstr  : non-linear constraint values of size m_nlcon (output)
 * m_ineq, Aineq, bineq : Aineq*x <= bineq constraint
 *             Aineq is an m_ineq-by-n matrix stored in row-major order (line by line)
 *             bineq is of size m_ineq
 * m_eq, Aeq, beq : Aeq*x = beq constraint
 *             Aeq is an m_eq-by-n matrix stored in row-major order (line by line)
 *             beq is of size m_eq
 * xl, xu    : x lower & upper bounds, of size n
 *
 * return    : see prima_rc enum for return codes
 */

PRIMAC_API
int prima_bobyqa(const prima_obj calfun, const int n, double x[], double *f,
                 const double xl[], const double xu[],
                 int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint);

PRIMAC_API
int prima_newuoa(const prima_obj calfun, const int n, double x[], double *f,
                 int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint);

PRIMAC_API
int prima_uobyqa(const prima_obj calfun, const int n, double x[], double *f,
                 int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int iprint);

PRIMAC_API
int prima_cobyla(const int m_nlcon, const prima_objcon calcfc, const int n, double x[], double *f,
                 double *cstrv, double nlconstr[],
                 const int m_ineq, const double Aineq[], const double bineq[],
                 const int m_eq, const double Aeq[], const double beq[],
                 const double xl[], const double xu[],
                 int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int iprint);

PRIMAC_API
int prima_lincoa(const prima_obj calfun, const int n, double x[], double *f,
                 double *cstrv,
                 const int m_ineq, const double Aineq[], const double bineq[],
                 const int m_eq, const double Aeq[], const double beq[],
                 const double xl[], const double xu[],
                 int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint);

#ifdef __cplusplus
}
#endif

#endif
