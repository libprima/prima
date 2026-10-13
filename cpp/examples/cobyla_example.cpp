// This is an example to illustrate the usage of the COBYLA solver.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include "prima/prima.hpp"

using namespace prima;
using namespace Eigen;

// Objective: f(x) = (x1-5)^2 + (x2-4)^2
double objective(const VectorXd& x) {
    return std::pow(x(0) - 5.0, 2) + std::pow(x(1) - 4.0, 2);
}

int main() {
    std::cout << "=== COBYLA Example ===" << std::endl;

    // Simple constrained optimization:
    // min  (x1-5)^2 + (x2-4)^2
    // s.t. x1^2 - 9 <= 0   (i.e., |x1| <= 3)

    VectorXd x0(2);
    x0 << 0.0, 0.0;

    // Nonlinear constraint: x1^2 - 9 <= 0
    auto cons_fun = [](const VectorXd& x) -> VectorXd {
        return VectorXd::Constant(1, x(0) * x(0) - 9.0);
    };
    VectorXd nlc_lb(1);
    nlc_lb << -std::numeric_limits<double>::infinity();
    VectorXd nlc_ub(1);
    nlc_ub << 0.0;
    NonlinearConstraint nlc(cons_fun, nlc_lb, nlc_ub);
    auto nlc_func = transform_constraint_function(nlc);

    MinimizeOptions opts;
    opts.quiet = false;  // Print progress
    opts.rhoend = 1e-6;
    opts.maxfun = 5000;

    auto result = minimize(objective, x0, "cobyla", nullptr, nullptr, &nlc_func, opts);

    std::cout << "\nResult:" << std::endl;
    std::cout << "  x = [" << result.x(0) << ", " << result.x(1) << "]" << std::endl;
    std::cout << "  f = " << result.fun << std::endl;
    std::cout << "  constraint x1^2-9 = " << (result.x(0) * result.x(0) - 9.0) << std::endl;
    std::cout << "  nfev = " << result.nfev << std::endl;

    // Check: x1 should be near 3, x2 near 4, f near 4
    if (std::abs(result.x(0) - 3.0) < 1e-2 &&
        std::abs(result.x(1) - 4.0) < 1e-2 &&
        std::abs(result.fun - 4.0) < 1e-2) {
        std::cout << "\nExample PASSED." << std::endl;
        return 0;
    } else {
        std::cerr << "\nExample FAILED." << std::endl;
        return 1;
    }
}
