// An example to illustrate the use of NEWUOA.


#include "prima/prima.h"
#include <math.h>
#include <stdio.h>


// Objective function
static void fun(const double x[], double *const f, const void *data)
{
    const double x1 = x[0];
    const double x2 = x[1];
    *f = pow(x1-5, 2) + pow(x2-4, 2);
    (void)data;
}


// Callback function
static void callback(const int n, const double x[], const double f, const int nf, const int tr, const double cstrv, const int m_nlcon, const double nlconstr[], bool *const terminate)
{
    (void)n;
    (void)cstrv;
    (void)m_nlcon;
    (void)nlconstr;
    printf("Best point so far: x = {%g, %g}, f = %g, nf = %d, tr = %d\n", x[0], x[1], f, nf, tr);
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
    problem.calfun = &fun;
    problem.x0 = x0;

    // Set up the options
    prima_options_t options;
    prima_init_options(&options);
    options.iprint = PRIMA_MSG_EXIT;
    options.rhoend = 1e-6;
    options.maxfun = 500*n;
    options.callback = &callback;

    // Call the solver
    prima_result_t result;
    const prima_rc_t rc = prima_minimize(PRIMA_NEWUOA, problem, options, &result);

    // Print the result
    printf("x* = {%g, %g}, rc = %d, msg = '%s', evals = %d\n", result.x[0], result.x[1], rc, result.message, result.nf);

    // Check the result
    int success = (fabs(result.x[0] - 5) > 2e-2 || fabs(result.x[1] - 4) > 2e-2);

    // Free the result
    prima_free_result(&result);

    return success;
}
