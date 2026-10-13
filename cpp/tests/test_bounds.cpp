#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <limits>

#include "prima/infos.hpp"
#include "prima/bounds.hpp"
#include "prima/project.hpp"
#include "prima/linear_constraints.hpp"
#include "prima/nonlinear_constraints.hpp"
#include "prima/linalg.hpp"
#include "prima/prima.hpp"
#include "test_common.hpp"

using namespace prima;
using namespace Eigen;
using prima::test::TestRunner;

void test_process_bounds(TestRunner& t) {
    auto [lb, ub] = process_bounds(nullptr, 3);
    t.check(lb.size() == 3, "process_bounds(nullptr, 3) size");
    t.check(std::isinf(-lb(0)) && lb(0) < 0, "process_bounds nullptr lb = -inf");
    t.check(std::isinf(ub(0)) && ub(0) > 0, "process_bounds nullptr ub = inf");

    VectorXd lb_in(2); lb_in << 1.0, 2.0;
    VectorXd ub_in(2); ub_in << 3.0, 4.0;
    Bounds b(lb_in, ub_in);
    auto [lb2, ub2] = process_bounds(&b, 3);
    t.check(lb2.size() == 3, "process_bounds bounds size");
    t.check(std::abs(lb2(0) - 1.0) < 1e-15, "process_bounds lb[0]");
    t.check(std::abs(lb2(1) - 2.0) < 1e-15, "process_bounds lb[1]");
    t.check(std::isinf(-lb2(2)) && lb2(2) < 0, "process_bounds lb[2] padded -inf");
    t.check(std::abs(ub2(0) - 3.0) < 1e-15, "process_bounds ub[0]");
    t.check(std::abs(ub2(1) - 4.0) < 1e-15, "process_bounds ub[1]");
    t.check(std::isinf(ub2(2)) && ub2(2) > 0, "process_bounds ub[2] padded inf");
}

void test_eliminate_fixed_bounds(TestRunner& t) {
    auto obj = [](const VectorXd& x) -> double { return x.squaredNorm(); };

    VectorXd lb(5), ub(5);
    lb << -1, -std::numeric_limits<double>::infinity(), 1,
           -std::numeric_limits<double>::infinity(), -0.5;
    ub << -0.5, -0.5, std::numeric_limits<double>::infinity(),
           std::numeric_limits<double>::infinity(), -0.5;
    Bounds bounds(lb, ub);
    VectorXd x0(5); x0 << 1, 2, 3, 4, 5;

    MinimizeOptions opts;
    opts.quiet = true;
    auto result = minimize(obj, x0, "cobyla", &bounds, nullptr, nullptr, opts);
    t.check_close(result.x(0), -0.5, 1e-3, "x[0] == -0.5");
    t.check_close(result.x(1), -0.5, 1e-3, "x[1] == -0.5");
    t.check_close(result.x(2), 1.0, 1e-3, "x[2] == 1");
    t.check_close(result.x(3), 0.0, 5e-3, "x[3] == 0");
    t.check_close(result.x(4), -0.5, 1e-3, "x[4] == -0.5");
    t.check_close(result.fun, 1.75, 1e-3, "f == 1.75");
}

void test_eliminate_fixed_bounds_with_linear_constraints(TestRunner& t) {
    auto obj = [](const VectorXd& x) -> double { return x.squaredNorm(); };

    VectorXd lb(3), ub(3);
    lb << -1, -std::numeric_limits<double>::infinity(),
           -std::numeric_limits<double>::infinity();
    ub << -1, std::numeric_limits<double>::infinity(),
           std::numeric_limits<double>::infinity();
    Bounds bounds(lb, ub);
    MatrixXd A(1, 3); A << 1, 1, 1;
    LinearConstraint lc(A, VectorXd::Constant(1, 9), VectorXd::Constant(1, 15));
    VectorXd x0(3); x0 << 1, 1, 1;

    MinimizeOptions opts;
    opts.quiet = true;
    auto result = minimize(obj, x0, "cobyla", &bounds, &lc, nullptr, opts);
    t.check_close(result.x(0), -1.0, 1e-6, "x[0] == -1");
    t.check_close(result.x(1), 5.0, 5e-5, "x[1] == 5");
    t.check_close(result.x(2), 5.0, 5e-5, "x[2] == 5");
    t.check_close(result.fun, 51.0, 1e-6, "f == 51");
}

void test_eliminate_fixed_bounds_with_nonlinear_constraints(TestRunner& t) {
    auto obj = [](const VectorXd& x) -> double { return x.squaredNorm(); };

    VectorXd lb(3), ub(3);
    lb << -1, -std::numeric_limits<double>::infinity(),
           -std::numeric_limits<double>::infinity();
    ub << -1, std::numeric_limits<double>::infinity(),
           std::numeric_limits<double>::infinity();
    Bounds bounds(lb, ub);

    int n = 3;
    auto nlc_fun = [n](const VectorXd& x) -> VectorXd {
        return VectorXd::Constant(1, x(n - 1) * x(n - 1));
    };
    NonlinearConstraint nlc(nlc_fun, VectorXd::Constant(1, 9), VectorXd::Constant(1, 15));
    auto nlc_func = transform_constraint_function(nlc);
    VectorXd x0(3); x0 << 1, 1, 1;

    MinimizeOptions opts;
    opts.quiet = true;
    auto result = minimize(obj, x0, "cobyla", &bounds, nullptr, &nlc_func, opts);
    t.check_close(result.x(0), -1.0, 1e-6, "x[0] == -1");
    t.check_close(result.x(1), 0.0, 1e-6, "x[1] == 0");
    t.check_close(result.x(2), 3.0, 1e-6, "x[2] == 3");
    t.check_close(result.fun, 10.0, 1e-6, "f == 10");
}

void test_all_fixed(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    Bounds bounds(Vector2d(1.0, 2.0), Vector2d(1.0, 2.0));
    auto obj = [](const VectorXd& x) -> double {
        return std::pow(x(0) - 1.0, 2) + std::pow(x(1) - 2.0, 2);
    };

    MinimizeOptions opts;
    opts.quiet = true;
    auto result = minimize(obj, x0, "cobyla", &bounds, nullptr, nullptr, opts);
    t.check(result.success, "all fixed success");
    t.check(result.status == FIXED_SUCCESS, "status FIXED_SUCCESS");
    t.check_close(result.x(0), 1.0, 1e-15, "x[0] == 1");
    t.check_close(result.x(1), 2.0, 1e-15, "x[1] == 2");
    t.check(result.nfev == 1, "nfev == 1");
}

int main() {
    std::cout << "=== test_bounds: Bounds and Fixed Variables ===" << std::endl;
    TestRunner t;

    test_process_bounds(t);
    test_eliminate_fixed_bounds(t);
    test_eliminate_fixed_bounds_with_linear_constraints(t);
    test_eliminate_fixed_bounds_with_nonlinear_constraints(t);
    test_all_fixed(t);

    return t.done();
}
