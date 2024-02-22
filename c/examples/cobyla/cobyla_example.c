// An example to illustrate the use of COBYLA.


#include "prima/prima.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>


#define M_NLCON 1
#define PROVIDE_INITIAL_F_AND_NLCONSTR 1


// Objective & constraint function
static void fun(const double x[], double *f, double constr[], const void *data)
{
    const double x1 = x[0];
    const double x2 = x[1];
    *f = pow(x1-5, 2) + pow(x2-4, 2);
    constr[0] = pow(x1, 2) - 9;  // Constraint: x1^2 - 9 <= 0
    (void)data;
}


// Callback function
static void callback(const int n, const double x[], const double f, const int nf, const int tr, const double cstrv, const int m_nlcon, const double nlconstr[], bool *terminate)
{
    (void)n;
    (void)m_nlcon;
    printf("Best point so far: x = {%g, %g}, f = %g, cstrv = %g, nlconstr = {%g}, nf = %d, tr = %d\n", x[0], x[1], f, cstrv, nlconstr[0], nf, tr);
    *terminate = 0;
}


// Main function
int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const int n = 2;
    double x0[2] = {0.0, 0.0};

    // Set up the problem
    prima_problem_t problem;
    prima_init_problem(&problem, n);
    problem.x0 = x0;
    problem.m_nlcon = M_NLCON;
    problem.calcfc = &fun;
    // Provide the initial values of the objective function and the nonlinear constraints.
    // This is OPTIONAL, and end users should NOT do it in general. Here, we do it for testing.
#if PROVIDE_INITIAL_F_AND_NLCONSTR
    double nlconstr0[M_NLCON] = {0};
    fun(x0, &(problem.f0), nlconstr0, NULL);
    problem.nlconstr0 = nlconstr0;
#endif

    // Set up the options
    prima_options_t options;
    prima_init_options(&options);
    options.iprint = PRIMA_MSG_EXIT;
    options.rhoend = 1e-6;
    options.maxfun = 500*n;
    options.callback = &callback;

    // Call the solver
    prima_result_t result;
    const int rc = prima_minimize(PRIMA_COBYLA, &problem, &options, &result);

    // Print the result
    printf("x* = {%g, %g}, f* = %g, cstrv = %g, nlconstr = {%g}, rc = %d, msg = '%s', evals = %d\n", result.x[0], result.x[1], result.f, result.cstrv, result.nlconstr[0], rc, result.message, result.nf);

    // Check the result
    int success = (fabs(result.x[0] - 3) > 2e-2 || fabs(result.x[1] - 4) > 2e-2);

    // Free the result
    prima_free_result(&result);

    return success;
}
