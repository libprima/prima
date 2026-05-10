// An example to illustrate the use of BOBYQA.


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
    problem.x0 = x0;
    problem.calfun = &fun;
    // We define an upper bound that will be active in order to demonstrate the usage.
    double xl[] = {-1.0, -1.0};
    double xu[] = {4.5, 4.5};
    problem.xl = xl;
    problem.xu = xu;

    // Set up the options
    prima_options_t options;
    prima_init_options(&options);
    options.iprint = PRIMA_MSG_EXIT;
    options.rhoend = 1e-6;
    options.maxfun = 500*n;
    options.callback = &callback;

    // Call the solver
    prima_result_t result;
    const prima_rc_t rc = prima_minimize(PRIMA_BOBYQA, problem, options, &result);

    // Print the result
    printf("x* = {%g, %g}, rc = %d, msg = '%s', evals = %d\n", result.x[0], result.x[1], rc, result.message, result.nf);

    // Check the result
    int success = (fabs(result.x[0] - 4.5) > 2e-2 || fabs(result.x[1] - 4) > 2e-2);

    // Free the result
    prima_free_result(&result);

    // for the second part of the example, we will change the initial point close to the upper bound
    // and set honour_x0 to true
    double new_x0[2] = {4.4, 4.4};

    // Set up the problem
    prima_problem_t honour_problem;
    prima_init_problem(&honour_problem, n);
    honour_problem.x0 = new_x0;
    honour_problem.calfun = &fun;
    // We define an upper bound that will be active in order to demonstrate the usage.
    honour_problem.xl = xl;
    honour_problem.xu = xu;

    prima_options_t honour_options;
    prima_init_options(&honour_options);
    honour_options.iprint = PRIMA_MSG_EXIT;
    honour_options.rhoend = 1e-6;
    honour_options.maxfun = 500*n;
    honour_options.honour_x0 = true; // change this to false to see the warning of resetting x0
    honour_options.callback = &callback;
    // Call the solver again with honour_x0 = true
    prima_result_t honour_result;
    const prima_rc_t honour_rc = prima_minimize(PRIMA_BOBYQA, honour_problem, honour_options, &honour_result);
    
    // Print the result
    printf("x* = {%g, %g}, rc = %d, msg = '%s', evals = %d\n", honour_result.x[0], honour_result.x[1], honour_rc, honour_result.message, honour_result.nf);

    // Check the result
    int honour_success = (fabs(honour_result.x[0] - 4.5) > 2e-2 || fabs(honour_result.x[1] - 4) > 2e-2);

    return honour_success;
}
