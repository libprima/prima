// Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).

#ifndef PRIMA_H
#define PRIMA_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#define PRIMAC_HELPER_DLL_IMPORT __declspec(dllimport)
#define PRIMAC_HELPER_DLL_EXPORT __declspec(dllexport)
#else
#define PRIMAC_HELPER_DLL_IMPORT
#define PRIMAC_HELPER_DLL_EXPORT
#endif

#ifdef PRIMAC_STATIC
#define PRIMAC_API
#elif defined primac_EXPORTS
#define PRIMAC_API PRIMAC_HELPER_DLL_EXPORT
#else
#define PRIMAC_API PRIMAC_HELPER_DLL_IMPORT
#endif

// Verbosity level
typedef enum {
  PRIMA_MSG_NONE = 0,  // Do not print any message
  PRIMA_MSG_EXIT = 1,  // Print a message at exit
  PRIMA_MSG_RHO = 2,  // Print a message when rho changes
  PRIMA_MSG_FEVL = 3,  // Print a message when the object/constraint functions get evaluated
} prima_message_t;


// Possible return values
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
  PRIMA_PROBLEM_SOLVER_MISMATCH_NONLINEAR_CONSTRAINTS = 115,
  PRIMA_PROBLEM_SOLVER_MISMATCH_LINEAR_CONSTRAINTS = 116,
  PRIMA_PROBLEM_SOLVER_MISMATCH_BOUNDS = 117,
} prima_rc_t;


// Function to get the message string corresponding to a return code
PRIMAC_API
const char *prima_get_rc_string(const prima_rc_t rc);


/*
 * Objective function required by UOBYQA, NEWUOA, BOBYQA, and LINCOA
 * x     : on input, the vector of variables (should not be modified)
 * f     : on output, the value of the function
 *         a NaN value can be passed to signal an evaluation error
 * data  : user data
*/
typedef void (*prima_obj_t)(const double x[], double *f, const void *data);


/*
 * Objective & constraint function required by COBYLA
 * x     : on input, the vector of variables (should not be modified)
 * f     : on output, the value of the function
 *         a NaN value can be passed to signal an evaluation error
 * constr : on output, the value of the constraints (of size m_nlcon),
 *          with the constraints being constr <= 0
 *          NaN values can be passed to signal evaluation errors
 * data  : user data
*/
typedef void (*prima_objcon_t)(const double x[], double *f, double constr[], const void *data);

/*
 * An optional callback function to report algorithm progress
 * n     : number of variables
 * x     : the current best point
 * f     : the function value of the current best point
 * nf    : number of function evaluations
 * tr    : number of trust-region iterations
 * cstrv : the constraint violation of the current best point (LINCOA and COBYLA only)
 * m_nlcon : number of nonlinear constraints (COBYLA only)
 * nlconstr : nonlinear constraint values of the current best point (COBYLA only)
 * terminate : a boolean to ask for termination
*/
typedef void (*prima_callback_t)(const int n, const double x[], const double f, int nf, int tr,
                               const double cstrv, int m_nlcon, const double nlconstr[], bool *terminate);

// Structure to hold the options
typedef struct {

  // rhobeg: a reasonable initial change to the variables; by default this is NaN, which will be
  // interpreted in Fortran as not present, in which case a default value will be used
  double rhobeg;

  // rhoend: required accuracy for the variables; by default this is NaN, which will be interpreted
  // in Fortran as not present, in which case a default value will be used
  double rhoend;

  // maxfun: maximum number of function evaluations; by default this is 0, which will be interpreted
  // in Fortran as not present, in which case a default value will be used
  int maxfun;

  // iprint: verbosity level, see the prima_message_t enum; default = PRIMA_MSG_NONE
  int iprint;

  // ftarget: target function value; solver stops when f <= ftarget for a feasible point; default: -Inf
  double ftarget;

  // npt: number of points in the interpolation set for NEWUOA, BOBYQA, and LINCOA; it should
  // satisfy n+2 <= npt <= (n+1)(n+2)/2; by default this is 0, which will be interpreted by Fortran
  // as not present, in which case a default value based on the algorithm will be used;
  // it is not used by UOBYQA or COBYLA.
  int npt;

  // user data, will be passed through the objective function callback; default: NULL
  void *data;

  // callback function to report algorithm progress; default: NULL
  prima_callback_t callback;

} prima_options_t;


// Function to initialize the options
PRIMAC_API
int prima_init_options(prima_options_t *options);


// Structure to hold the problem
typedef struct {

  // dimension of the problem
  int n;

  // objective function to minimize (not COBYLA)
  prima_obj_t calfun;

  // objective function to minimize with constraints (COBYLA only)
  prima_objcon_t calcfc;

  // starting point
  double *x0;

  // bound constraints: xl <= x <= xu
  double *xl;
  double *xu;

  // linear inequality constraints: Aineq*x <= bineq
  // Aineq is an m_ineq-by-n matrix stored in row-major order (row by row)
  // bineq is of size m_ineq
  int m_ineq;
  double *Aineq;
  double *bineq;

  // linear equality constraints: Aeq*x = beq
  // Aeq is an m_eq-by-n matrix stored in row-major order (row by row)
  // beq is of size m_eq
  int m_eq;
  double *Aeq;
  double *beq;

  // number of nonlinear constraints for COBYLA; default: 0
  int m_nlcon;

  // initial objective function value and constraint values (COBYLA only)
  // They are OPTIONAL and should ONLY be used when interfacing with MATLAB/Python/Julia/R,
  // where we need to evaluate the constraints at the initial point to get m_nlcon.
  // C end users should leave them as default, i.e., f0 = NAN and nlconstr0 = NULL.
  double f0;
  double *nlconstr0;

} prima_problem_t;


// Function to initialize the problem
PRIMAC_API
int prima_init_problem(prima_problem_t *problem, int n);


// Structure to hold the result
typedef struct {

  // returned point
  double *x;

  // objective function value at the returned point
  double f;

  // number of objective function evaluations
  int nf;

  // constraint violation at the returned point (COBYLA and LINCOA only)
  double cstrv;

  // nonlinear constraint values at the returned point, of size m_nlcon (COBYLA only)
  double *nlconstr;

  // exit code
  int status;

  // exit message
  const char *message;

} prima_result_t;


// Function to free the result
PRIMAC_API
int prima_free_result(prima_result_t *result);


// Algorithm
typedef enum
{
  PRIMA_BOBYQA,
  PRIMA_COBYLA,
  PRIMA_LINCOA,
  PRIMA_NEWUOA,
  PRIMA_UOBYQA
} prima_algorithm_t;


/*
 * The function that does the minimization using a PRIMA solver
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
