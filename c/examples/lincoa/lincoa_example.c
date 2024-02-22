// An example to illustrate the use of LINCOA.


#include "prima/prima.h"
#include <math.h>
#include <stdio.h>


// Objective function
static void fun(const double x[], double *f, const void *data)
{
    const double x1 = x[0];
    const double x2 = x[1];
    *f = pow(x1-5, 2) + pow(x2-4, 2);
    (void)data;
}


// Callback function
static void callback(int n, const double x[], double f, int nf, int tr, double cstrv, int m_nlcon, const double nlconstr[], bool *terminate)
{
    (void)n;
    (void)m_nlcon;
    (void)nlconstr;
    printf("Best point so far: x = {%g, %g}, f = %g, cstrv = %g, nf = %d, tr = %d\n", x[0], x[1], f, cstrv, nf, tr);
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
    // We define linear constraints that will be active in order to demonstrate the usage.
    // The constraint is x1 + x2 <= 8.
    problem.m_ineq = 1;
    double Aineq[1*2] = {1.0, 1.0};
    double bineq[1] = {8.0};
    problem.Aineq = Aineq;
    problem.bineq = bineq;

    // Set up the options
    prima_options_t options;
    prima_init_options(&options);
    options.iprint = PRIMA_MSG_EXIT;
    options.rhoend = 1e-6;
    options.maxfun = 500*n;
    options.callback = &callback;

    // Call the solver
    prima_result_t result;
    const int rc = prima_minimize(PRIMA_LINCOA, &problem, &options, &result);

    // Print the result
    printf("x* = {%g, %g}, f* = %g, cstrv = %g, rc = %d, msg = '%s', evals = %d\n", result.x[0], result.x[1], result.f, result.cstrv, rc, result.message, result.nf);

    // Check the result
    int success = (fabs(result.x[0] - 4.5) > 2e-2 || fabs(result.x[1] - 3.5) > 2e-2);

    // Free the result
    prima_free_result(&result);

    return success;
}
