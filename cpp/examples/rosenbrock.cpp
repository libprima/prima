// Illustration of how to use prima with COBYLA.
//
// Minimize the chained Rosenbrock function subject to various constraints.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <iostream>
#include <cmath>
#include <iomanip>
#include <Eigen/Core>
#include "prima/prima.hpp"

using namespace prima;
using namespace Eigen;

// Chained Rosenbrock function
double chrosen(const VectorXd& x) {
    int n = x.size();
    double f = 0;
    for (int i = 0; i < n - 1; ++i) {
        f += std::pow(1 - x(i), 2) + 4 * std::pow(x(i + 1) - x(i) * x(i), 2);
    }
    return f;
}

// Nonlinear inequality constraint: x(i)^2 >= x(i+1)
VectorXd nlc_ineq(const VectorXd& x) {
    int n = x.size();
    VectorXd c(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        c(i) = x(i) * x(i) - x(i + 1);
    }
    return c;
}

// Nonlinear equality constraint: ||x||^2 = 1
VectorXd nlc_eq(const VectorXd& x) {
    VectorXd c(1);
    c(0) = x.squaredNorm() - 1;
    return c;
}

int main() {
    std::cout << std::setprecision(4);
    std::cout << "Minimize the chained Rosenbrock function with three variables "
              << "subject to various constraints using COBYLA.\n" << std::endl;

    VectorXd x0(3);
    x0 << 0, 0, 0;

    // ---------------------------------------------------------------- //
    // 1. Nonlinear constraints
    //    ||x||_2^2 = 1, x(i)^2 >= x(i+1) >= 0.5*x(i) >= 0 for i = 1, 2
    // ---------------------------------------------------------------- //
    std::cout << "1. Nonlinear constraints --- ||x||_2^2 = 1, "
              << "x(i)^2 >= x(i+1) >= 0.5*x(i) >= 0 for i = 1, 2:\n" << std::endl;

    VectorXd lb(3); lb << 0, 0, 0;
    VectorXd ub(3); ub << std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity();
    Bounds bounds(lb, ub);

    // Linear constraints: 0.5*x(i) - x(i+1) <= 0
    MatrixXd A(2, 3);
    A << 0.5, -1,   0,
          0,    0.5, -1;
    LinearConstraint lin_con(A,
        VectorXd::Constant(2, -std::numeric_limits<double>::infinity()),
        VectorXd::Constant(2, 0.0));

    // Nonlinear constraints
    NonlinearConstraint nlc_ineq_obj(nlc_ineq,
        VectorXd::Constant(2, 0.0),
        VectorXd::Constant(2, std::numeric_limits<double>::infinity()));
    NonlinearConstraint nlc_eq_obj(nlc_eq,
        VectorXd::Constant(1, 0.0),
        VectorXd::Constant(1, 0.0));

    auto nlc_ineq_t = transform_constraint_function(nlc_ineq_obj);
    auto nlc_eq_t = transform_constraint_function(nlc_eq_obj);

    MinimizeOptions opts;
    opts.quiet = true;

    // COBYLA only handles one NonlinearConstraintFunction, so combine them
    NonlinearConstraintFunction combined_nlc = [nlc_ineq_t, nlc_eq_t](const VectorXd& x) -> VectorXd {
        VectorXd v1 = nlc_ineq_t(x);
        VectorXd v2 = nlc_eq_t(x);
        VectorXd r(v1.size() + v2.size());
        r.head(v1.size()) = v1;
        r.tail(v2.size()) = v2;
        return r;
    };

    auto result = minimize(chrosen, x0, "cobyla", &bounds, &lin_con, &combined_nlc, opts);
    std::cout << "  x = [" << result.x(0) << ", " << result.x(1) << ", " << result.x(2) << "]" << std::endl;
    std::cout << "  f = " << result.fun << "  nfev = " << result.nfev << std::endl;
    std::cout << "  ||x||^2 = " << result.x.squaredNorm() << std::endl;

    // ---------------------------------------------------------------- //
    // 2. Linear constraints
    //    sum(x) = 1, x(i+1) <= x(i) <= 1 for i = 1, 2
    // ---------------------------------------------------------------- //
    std::cout << "\n2. Linear constraints --- sum(x) = 1, x(i+1) <= x(i) <= 1 for i = 1, 2:\n" << std::endl;

    Bounds bounds2(
        VectorXd::Constant(3, -std::numeric_limits<double>::infinity()),
        VectorXd::Constant(3, 1.0));
    MatrixXd A2(3, 3);
    A2 << -1,  1,  0,
           0, -1,  1,
           1,  1,  1;
    LinearConstraint lin_con2(A2,
        Vector3d(-std::numeric_limits<double>::infinity(),
                 -std::numeric_limits<double>::infinity(), 1.0),
        Vector3d(0.0, 0.0, 1.0));

    auto result2 = minimize(chrosen, x0, "cobyla", &bounds2, &lin_con2, nullptr, opts);
    std::cout << "  x = [" << result2.x(0) << ", " << result2.x(1) << ", " << result2.x(2) << "]" << std::endl;
    std::cout << "  f = " << result2.fun << "  nfev = " << result2.nfev << std::endl;
    std::cout << "  sum(x) = " << result2.x.sum() << std::endl;

    // ---------------------------------------------------------------- //
    // 3. Bound constraints: -0.5 <= x(1) <= 0.5, 0 <= x(2) <= 0.25
    // ---------------------------------------------------------------- //
    std::cout << "\n3. Bound constraints --- -0.5 <= x(1) <= 0.5, 0 <= x(2) <= 0.25:\n" << std::endl;

    Bounds bounds3(
        Vector3d(-0.5, 0.0, -std::numeric_limits<double>::infinity()),
        Vector3d(0.5, 0.25, std::numeric_limits<double>::infinity()));

    auto result3 = minimize(chrosen, x0, "cobyla", &bounds3, nullptr, nullptr, opts);
    std::cout << "  x = [" << result3.x(0) << ", " << result3.x(1) << ", " << result3.x(2) << "]" << std::endl;
    std::cout << "  f = " << result3.fun << "  nfev = " << result3.nfev << std::endl;

    // ---------------------------------------------------------------- //
    // 4. No constraints
    // ---------------------------------------------------------------- //
    std::cout << "\n4. No constraints:\n" << std::endl;
    auto result4 = minimize(chrosen, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    std::cout << "  x = [" << result4.x(0) << ", " << result4.x(1) << ", " << result4.x(2) << "]" << std::endl;
    std::cout << "  f = " << result4.fun << "  nfev = " << result4.nfev << std::endl;

    return 0;
}
