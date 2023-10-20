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
  PRIMA_NULL_OPTIONS = 200,
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
 * data  : user-data
 * constr : on output, the value of the constraints (of size m_nlcon)
 *          NaN values can be passed to signal evaluation errors
 *          only for cobyla
*/
typedef void (*prima_obj)(const double x[], double *f, const void *data);
typedef void (*prima_objcon)(const double x[], double *f, double constr[], const void *data);


typedef struct {
  
  // a reasonable initial change to the variables (default=1.0)
  double rhobeg;

  // required accuracy for the variables (default=1e-6)
  double rhoend;

  // maximum number of function evaluations (default=100)
  int maxfun;

  // verbosity level, see the prima_message enum (default=PRIMA_MSG_NONE)
  int iprint;

  // target function value; optimization stops when f <= ftarget for a feasible point (default=-inf)
  double ftarget;

  // number of points in the interpolation set n+2<=npt<=(n+1)(n+2)/2 (default=2*n+1)
  int npt;

  // user-data, will be passed through the objective function callback
  void *data;

  // Aineq*x <= bineq constraint
  // Aineq is an m_ineq-by-n matrix stored in row-major order (line by line)
  // bineq is of size m_ineq
  int m_ineq;
  double *Aineq;
  double *bineq;

  // m_eq, Aeq, beq : Aeq*x = beq constraint
  // Aeq is an m_eq-by-n matrix stored in row-major order (line by line)
  // beq is of size m_eq
  int m_eq;
  double *Aeq;
  double *beq;

} prima_options;

/* Initialize option data */
PRIMAC_API
int prima_init_options(prima_options * options);

/*
 * calfun    : function to minimize (see prima_obj)
 * n         : number of variables (>=0)
 * x         : on input, initial estimate
 *             on output, the solution
 * f         : objective value (output)
 * nf        : number of objective function calls (output)
 * m_nlcon   : number of non-linear constraints (>=0)
 * calcfc    : function to minimize and constraints (see prima_objcon)
 * cstrv     : constraint violation (output)
 * nlconstr  : non-linear constraint values of size m_nlcon (output)
 * xl, xu    : x lower & upper bounds, of size n
 *
 * return    : see prima_rc enum for return codes
 */

PRIMAC_API
int prima_bobyqa(const prima_obj calfun, const int n, double x[], double *f,
                 const double xl[], const double xu[],
                 int *nf, const prima_options * options);

PRIMAC_API
int prima_newuoa(const prima_obj calfun, const int n, double x[], double *f,
                 int *nf, const prima_options * options);

PRIMAC_API
int prima_uobyqa(const prima_obj calfun, const int n, double x[], double *f,
                 int *nf, const prima_options * options);

PRIMAC_API
int prima_cobyla(const int m_nlcon, const prima_objcon calcfc, const int n, double x[], double *f,
                 double *cstrv, double nlconstr[],
                 const double xl[], const double xu[],
                 int *nf, const prima_options * options);

PRIMAC_API
int prima_lincoa(const prima_obj calfun, const int n, double x[], double *f,
                 double *cstrv,
                 const double xl[], const double xu[],
                 int *nf, const prima_options * options);

#ifdef __cplusplus
}
#endif

#endif
