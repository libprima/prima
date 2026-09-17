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


// Possible algorithms
typedef enum {
    PRIMA_UOBYQA,  // Unconstrained
    PRIMA_NEWUOA,  // Unconstrained
    PRIMA_BOBYQA,  // Bound-constrained
    PRIMA_LINCOA,  // Linearly constrained
    PRIMA_COBYLA,  // Nonlinearly constrained
} prima_algorithm_t;


// Verbosity level
typedef enum {
    PRIMA_MSG_NONE = 0,  // Do not print any message
    PRIMA_MSG_EXIT = 1,  // Print a message at exit
    PRIMA_MSG_RHO = 2,  // Print a message when rho changes
    PRIMA_MSG_FEVL = 3,  // Print a message when the object/constraint functions get evaluated
} prima_message_t;


// Possible return values
typedef enum {
    PRIMA_RC_DFT = 0,
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
    PRIMA_RESULT_INITIALIZED = 115,
} prima_rc_t;


/*
 * Objective function required by UOBYQA, NEWUOA, BOBYQA, and LINCOA
 * x     : on input, the vector of variables (should not be modified)
 * f     : on output, the value of the function
 *         a NaN value can be passed to signal an evaluation error
 * data  : user data
 */
typedef void (*prima_obj_t)(const double x[], double *const f, const void *data);


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
typedef void (*prima_objcon_t)(const double x[], double *const f, double constr[], const void *data);


/*
 * Callback function to report algorithm progress
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
typedef void (*prima_callback_t)(const int n, const double x[], const double f, const int nf, const int tr,
                                 const double cstrv, const int m_nlcon, const double nlconstr[], bool *const terminate);


// Structure to hold the problem
// In the following, "Default" refers to the value set by `prima_init_problem`.
typedef struct {

    // n: number of variables, n >= 1
    // Default: 0
    int n;

    // calfun: pointer to the objective function to minimize
    // Should not be NULL for UOBYQA, NEWUOA, BOBYQA, and LINCOA; should be NULL for COBYLA
    // Default: NULL
    prima_obj_t calfun;

    // calcfc: pointer to the objective & constraint function to minimize
    // Should not be NULL for COBYLA; should be NULL for UOBYQA, NEWUOA, BOBYQA, and LINCOA
    // Default: NULL
    prima_objcon_t calcfc;

    // starting point
    // Should not be NULL
    // Default: NULL
    double *x0;

    // Bound constraints: xl <= x <= xu
    // Should be NULL for UOBYQA and NEWUOA
    // Default: xl = NULL and xu = NULL
    double *xl;
    double *xu;

    // Linear inequality constraints: Aineq*x <= bineq
    // Aineq is an m_ineq-by-n matrix stored in row-major order (row by row)
    // bineq is of size m_ineq
    // Should be m_ineq = 0, Aineq = NULL, and bineq = NULL for UOBYQA, NEWUOA, and BOBYQA
    // Default: m_ineq = 0, Aineq = NULL, and bineq = NULL
    int m_ineq;
    double *Aineq;
    double *bineq;

    // Linear equality constraints: Aeq*x = beq
    // Aeq is an m_eq-by-n matrix stored in row-major order (row by row)
    // beq is of size m_eq
    // Should be m_eq = 0, Aeq = NULL, and beq = NULL for UOBYQA, NEWUOA, and BOBYQA
    // Default: m_eq = 0, Aeq = NULL, and beq = NULL
    int m_eq;
    double *Aeq;
    double *beq;

    // m_nlcon: number of nonlinear constraints defined by calcfc
    // Should be 0 for UOBYQA, NEWUOA, BOBYQA, and LINCOA
    // Default: 0
    int m_nlcon;

    // f0, nlconstr0: initial objective function value and constraint values (COBYLA only)
    // It should ONLY be used when interfacing with high-level languages such as MATLAB/Python/
    // Julia/R. In these languages, instead of asking the user to provide m_nlcon, we should
    // evaluate the constraints at the initial point to get m_nlcon; we evaluate the objective
    // function at the initial point as well to keep the objective and constraint evaluations
    // synchronized; after this, we pass the initial object and constraint values to the solver
    // via f0 and nlconstr0 to avoid re-evaluating them, as the evaluations are expensive.
    // C end users should leave f0 and nlconstr0 as the default set by prima_init_problem and
    // set m_nlcon to the number of nonlinear constraints.
    // Default: f0 = NAN and nlconstr0 = NULL
    double f0;
    double *nlconstr0;

} prima_problem_t;


// Function to initialize the problem
PRIMAC_API
prima_rc_t prima_init_problem(prima_problem_t *const problem, const int n);


// Structure to hold the options
// In the following, "Default" refers to the value set by `prima_init_options`.
typedef struct {

    // rhobeg: a reasonable initial change to the variables
    // Default: NaN, which will be interpreted in Fortran as not present, in which case a default
    // value will be used
    double rhobeg;

    // rhoend: required accuracy for the variables
    // Default: NaN, which will be interpreted in Fortran as not present, in which case a default
    // value will be used
    double rhoend;

    // maxfun: maximum number of function evaluations
    // Default: 0, which will be interpreted in Fortran as not present, in which case a default value
    // will be used
    int maxfun;

    // iprint: verbosity level (see the prima_message_t enum)
    // Default: PRIMA_MSG_NONE, which means that no message will be printed
    int iprint;

    // ftarget: target function value; solver stops when f <= ftarget for a feasible point
    // Default: -Inf
    double ftarget;

    // npt: number of points in the interpolation set for NEWUOA, BOBYQA, and LINCOA
    // Should satisfy n+2 <= npt <= (n+1)(n+2)/2
    // It will be ignored by UOBYQA or COBYLA if provided
    // Default: 0, which will be interpreted by Fortran as not present, in which case a default value
    // based on the algorithm that will be used
    int npt;

    // ctol: tolerance for the constraint violation (COBYLA and LINCOA only)
    // ctol is the tolerance of constraint violation. x is considered feasible if cstrv(x) <= ctol.
    // N.B.: 1. ctol is absolute, not relative.
    // 2. ctol is used for choosing the returned x. It does not affect the iterations of the algorithm.
    // Default: NaN, which will be interpreted in Fortran as not present, in which case a default value
    // of sqrt(machine epsilon) will be used.
    double ctol;

    // data: user data, will be passed through the objective function
    // Default: NULL
    void *data;

    // callback: pointer to the callback function to report algorithm progress
    // Default: NULL, which means that no callback will be called
    prima_callback_t callback;

} prima_options_t;


// Function to initialize the options
PRIMAC_API
prima_rc_t prima_init_options(prima_options_t *const options);


// Functions implemented in Fortran (*_c.f90)
int cobyla_c(const int m_nlcon, const prima_objcon_t calcfc, const void *data, const int n, double x[], double *const f, double *const cstrv, double nlconstr[],
             const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[],
             const double xl[], const double xu[],
             const double f0, const double nlconstr0[],
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int iprint, const double ctol,
             const prima_callback_t callback, int *const info);

int bobyqa_c(prima_obj_t calfun, const void *data, const int n, double x[], double *const f, const double xl[], const double xu[],
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, const prima_callback_t callback, int *const info);

int newuoa_c(prima_obj_t calfun, const void *data, const int n, double x[], double *const f,
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, const prima_callback_t callback, int *const info);

int uobyqa_c(prima_obj_t calfun, const void *data, const int n, double x[], double *const f,
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int iprint, const prima_callback_t callback, int *const info);

int lincoa_c(prima_obj_t calfun, const void *data, const int n, double x[], double *const f,
             double *const cstrv, const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[], const double xl[], const double xu[],
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, const double ctol,
             const prima_callback_t callback, int *const info);

// Structure to hold the result
// prima_minimize will allocate the memory for x and nlconstr (if needed),
// and as such the user is expected to free the memory using prima_free_result
// once they are done using the contents of the result (or have saved the contents).
typedef struct {

    // x: returned point
    double *x;

    // f: objective function value at the returned point
    double f;

    // cstrv: constraint violation at the returned point (COBYLA and LINCOA only)
    double cstrv;

    // nlconstr: nonlinear constraint values at the returned point, of size m_nlcon (COBYLA only)
    double *nlconstr;

    // nf: number of function evaluations
    int nf;

    // status: return code
    prima_rc_t status;

    // success: whether the solver returned normally or ran into abnormal conditions
    bool success;

    // message: exit message
    const char *message;

} prima_result_t;


// Function to free the result
PRIMAC_API
prima_rc_t prima_free_result(prima_result_t *const result);


/*
 * The function that does the minimization using a PRIMA solver.
 * After using (or saving) the result, the user is expected to call
 * prima_free_result, regardless of the return value of prima_minimize,
 * in order to avoid memory leaks.
 * algorithm : optimization algorithm (see prima_algorithm_t)
 * problem   : optimization problem (see prima_problem_t)
 * options   : optimization options (see prima_options_t)
 * result    : optimization result (see prima_result_t)
 * return    : see prima_rc_t enum for return codes
 */
PRIMAC_API
prima_rc_t prima_minimize(const prima_algorithm_t algorithm, const prima_problem_t problem, const prima_options_t options, prima_result_t *const result);


#ifdef __cplusplus
}
#endif

#endif
