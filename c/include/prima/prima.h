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
  PRIMA_NULL_OPTIONS = 110,
  PRIMA_NULL_RESULT = 111,
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

  // maximum number of function evaluations (default=-1 interpreted as 500*n)
  int maxfun;

  // verbosity level, see the prima_message enum (default=PRIMA_MSG_NONE)
  int iprint;

  // target function value; optimization stops when f <= ftarget for a feasible point (default=-inf)
  double ftarget;

  // number of points in the interpolation set n+2<=npt<=(n+1)(n+2)/2 (default=-1 interpreted as 2*n+1)
  // ignored for uobyqa & cobyla
  int npt;

  // user-data, will be passed through the objective function callback
  void *data;

  // bound constraints, ignored for newuoa & uobyqa
  double *xl;
  double *xu;

  // whether prima had to allocate xl/xu (private, do not use)
  int _allocated_xl;
  int _allocated_xu;

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

  // number of non-linear constraints for cobyla (>=0), cobyla-only, default=0
  int m_nlcon;

  // should be set to the objective function value and constraints values of the starting X, cobyla-only
  double f0;
  double *nlconstr0;

  // whether prima had to allocate nlconstr0 (private, do not use)
  int _allocated_nlconstr0;
  
} prima_options;


/* Initialize/free option data */
PRIMAC_API
int prima_init_options(prima_options * options);

PRIMAC_API
int prima_free_options(prima_options * opt);

typedef struct {
  // objective value
  double f;

  // number of objective function calls
  int nf;

  // constraint violation (cobyla & lincoa)
  double cstrv;

  // non-linear constraint values, of size m_nlcon (cobyla only)
  double *nlconstr;

  // size of nlconstr (private, do not use)
  int _m_nlcon;

} prima_result;


/* Free result after optimization */
PRIMAC_API
int prima_free_result(prima_result * result);

/*
 * calfun    : function to minimize (see prima_obj)
 * calcfc    : function to minimize and constraints (see prima_objcon)
 * n         : number of variables (>=0)
 * x         : on input, initial estimate
 *             on output, the solution
 * options   : optimization options (see prima_options)
 * result    : optimization result (see prima_result)
 * return    : see prima_rc enum for return codes
 */

PRIMAC_API
int prima_bobyqa(const prima_obj calfun, const int n, double x[], prima_options *options, prima_result *result);

PRIMAC_API
int prima_newuoa(const prima_obj calfun, const int n, double x[], prima_options *options, prima_result *result);

PRIMAC_API
int prima_uobyqa(const prima_obj calfun, const int n, double x[], prima_options *options, prima_result *result);

PRIMAC_API
int prima_cobyla(const prima_objcon calcfc, const int n, double x[], prima_options *options, prima_result *result);

PRIMAC_API
int prima_lincoa(const prima_obj calfun, const int n, double x[], prima_options *options, prima_result *result);

#ifdef __cplusplus
}
#endif

#endif
