// Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).

#ifndef PRIMA_H
#define PRIMA_H

#include <stdbool.h>

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
} prima_message_t;

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
  PRIMA_CALLBACK_TERMINATE = 30,
  PRIMA_INVALID_INPUT = 100,
  PRIMA_ASSERTION_FAILS = 101,
  PRIMA_VALIDATION_FAILS = 102,
  PRIMA_MEMORY_ALLOCATION_FAILS = 103,
  PRIMA_NULL_OPTIONS = 110,
  PRIMA_NULL_PROBLEM = 111,
  PRIMA_NULL_X0 = 112,
  PRIMA_NULL_RESULT = 113,
  PRIMA_NULL_FUNCTION = 114,
} prima_rc_t;

/*
 * Return code string
 */
PRIMAC_API
const char *prima_get_rc_string(const prima_rc_t rc);

/*
 * The objective function as required by solvers
 *
 * x     : on input, then vector of variables (should not be modified)
 * f     : on output, the value of the function
 *         a NaN value can be passed to signal an evaluation error
 * data  : user data
 * constr : on output, the value of the constraints (of size m_nlcon)
 *          NaN values can be passed to signal evaluation errors
 *          only for cobyla
*/
typedef void (*prima_obj_t)(const double x[], double *f, const void *data);
typedef void (*prima_objcon_t)(const double x[], double *f, double constr[], const void *data);

/* An optional callback function to report algorithm progress
 *
 * n     : number of variables
 * x     : the current best point
 * f     : the function value of the best point
 * nf    : number of objective function calls
 * tr    : iteration number
 * cstrv : the constraint value verified by the current best point
 * m_nlcon : number of non-linear constraints (cobyla only)
 * nlconstr : non-linear constraints values verified by the current best point (cobyla only)
 * terminate : a boolean to ask from early optimization exit
*/
typedef void (*prima_callback_t)(const int n, const double x[], const double f, int nf, int tr,
                               const double cstrv, int m_nlcon, const double nlconstr[], bool *terminate);


typedef struct {
  
  // a reasonable initial change to the variables (default=1.0)
  double rhobeg;

  // required accuracy for the variables (default=1e-6)
  double rhoend;

  // maximum number of function evaluations (default=-1 interpreted as 500*n)
  int maxfun;

  // verbosity level, see the prima_message_t enum (default=PRIMA_MSG_NONE)
  int iprint;

  // target function value; optimization stops when f <= ftarget for a feasible point (default=-inf)
  double ftarget;

  // number of points in the interpolation set n+2<=npt<=(n+1)(n+2)/2 (default=-1 interpreted as 2*n+1)
  // ignored for uobyqa & cobyla
  int npt;

  // user data, will be passed through the objective function callback
  void *data;

  // callback function to report algorithm progress (default=NULL)
  prima_callback_t callback;

} prima_options_t;

/* Initialize problem */
PRIMAC_API
int prima_init_options(prima_options_t *options);

typedef struct {

  // dimension of the problem
  int n;

  // objective function to minimize (not cobyla)
  prima_obj_t calfun;

  // objective function to minimize with constraints (cobyla)
  prima_objcon_t calcfc;
  
  // starting point
  double *x0;

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
  
} prima_problem_t;


/* Initialize/free problem */
PRIMAC_API
int prima_init_problem(prima_problem_t *problem, int n);

PRIMAC_API
int prima_free_problem(prima_problem_t *problem);


typedef struct {

  // final point
  double *x;

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

  // exit code
  int status;
  
  // error message
  const char *message;

} prima_result_t;


/* Free result after optimization */
PRIMAC_API
int prima_free_result(prima_result_t * result);

/*
 * Algorithm
 */
typedef enum
{
  PRIMA_BOBYQA,
  PRIMA_COBYLA,
  PRIMA_LINCOA,
  PRIMA_NEWUOA,
  PRIMA_UOBYQA
} prima_algorithm_t;


/*
 * algorithm : optimization algorithm (see prima_algorithm)
 * problem   : optimization problem (see prima_problem)
 * options   : optimization options (see prima_options)
 * result    : optimization result (see prima_result)
 * return    : see prima_rc_t enum for return codes
 */

PRIMAC_API
int prima_minimize(const prima_algorithm_t algorithm, prima_problem_t *problem, prima_options_t *options, prima_result_t *result);

#ifdef __cplusplus
}
#endif

#endif
