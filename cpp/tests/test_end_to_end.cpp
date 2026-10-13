#include <Eigen/Core>
#include <cmath>
#include <iostream>

#include "prima/infos.hpp"
#include "prima/bounds.hpp"
#include "prima/linear_constraints.hpp"
#include "prima/nonlinear_constraints.hpp"
#include "prima/prima.hpp"
#include "test_common.hpp"

using namespace prima;
using namespace Eigen;
using prima::test::TestRunner;

static double obj(const VectorXd& x) {
    return std::pow(x(0) - 1.0, 2) + std::pow(x(1) - 2.5, 2);
}

void test_no_constraints(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    MinimizeOptions opts;
    opts.quiet = true;
    opts.rhoend = 1e-4;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "no constraints success");
    t.check_close(result.x(0), 1.0, 1e-3, "x[0] close to 1");
    t.check_close(result.x(1), 2.5, 1e-3, "x[1] close to 2.5");
}

void test_bounds(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    Bounds bounds(Vector2d(-5.0, 10.0), Vector2d(5.0, 10.0));
    MinimizeOptions opts;
    opts.quiet = true;
    opts.rhoend = 1e-4;
    auto result = minimize(obj, x0, "cobyla", &bounds, nullptr, nullptr, opts);
    t.check(result.success, "bounds success");
    t.check_close(result.x(0), 1.0, 1e-3, "x[0] close to 1");
    t.check_close(result.x(1), 10.0, 1e-3, "x[1] close to 10");
}

void test_linear_constraints(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    // x0 + x1 = 5
    MatrixXd A(1, 2); A << 1.0, 1.0;
    LinearConstraint lc(A, VectorXd::Constant(1, 5.0), VectorXd::Constant(1, 5.0));
    MinimizeOptions opts;
    opts.quiet = true;
    opts.rhoend = 1e-4;
    auto result = minimize(obj, x0, "cobyla", nullptr, &lc, nullptr, opts);
    t.check(result.success, "linear constraints success");
    t.check_close(result.x(0), 1.75, 1.0, "x[0] close to 1.75");
    t.check_close(result.x(1), 3.25, 1.0, "x[1] close to 3.25");
    t.check(std::abs(result.x(0) + result.x(1) - 5.0) < 1e-3, "eq satisfied");
}

void test_nonlinear_constraint(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;

    auto cons_fun = [](const VectorXd& x) -> VectorXd {
        return VectorXd::Constant(1, x(0) * x(0) - 0.25);
    };
    NonlinearConstraint nlc(cons_fun, VectorXd::Constant(1, 0.0), VectorXd::Constant(1, 0.0));
    auto nlc_func = transform_constraint_function(nlc);

    MinimizeOptions opts;
    opts.quiet = true;
    opts.rhoend = 1e-4;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, &nlc_func, opts);
    t.check(result.success, "nonlinear constraint success");
    t.check_close(result.x(0), 0.5, 1e-3, "x[0] close to 0.5");
    t.check_close(result.x(1), 2.5, 0.5, "x[1] close to 2.5");
    t.check(std::abs(result.x(0) * result.x(0) - 0.25) < 1e-8, "cons satisfied");
}

void test_auto_detect(TestRunner& t) {
    VectorXd x0(2); x0 << 3.0, 3.0;
    MinimizeOptions opts;
    opts.quiet = true;
    opts.rhoend = 1e-4;
    opts.maxfun = 100;
    auto result = minimize(obj, x0, "", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "auto-detect success");
    t.check(result.method == "cobyla", "method == cobyla");
}

int main() {
    std::cout << "=== test_end_to_end: Solver Integration ===" << std::endl;
    TestRunner t;

    test_no_constraints(t);
    test_bounds(t);
    test_linear_constraints(t);
    test_nonlinear_constraint(t);
    test_auto_detect(t);

    return t.done();
}
