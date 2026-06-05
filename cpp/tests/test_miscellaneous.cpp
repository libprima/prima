#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <limits>

#include "prima/infos.hpp"
#include "prima/nonlinear_constraints.hpp"
#include "prima/prima.hpp"
#include "test_common.hpp"

using namespace prima;
using namespace Eigen;
using prima::test::TestRunner;

static double obj(const VectorXd& x) {
    return std::pow(x(0) - 1.0, 2) + std::pow(x(1) - 2.5, 2);
}

void test_callback_terminate(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    int call_count = 0;
    auto cb = [&](const VectorXd&, double, int, int, double, const VectorXd&) -> bool {
        ++call_count; return true;
    };
    MinimizeOptions opts;
    opts.quiet = true;
    opts.callback = cb;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(call_count >= 1, "callback was called");
    t.check(result.fun < 100, "non-catastrophic result");
}

void test_callback_no_terminate(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    int call_count = 0;
    auto cb = [&](const VectorXd&, double, int, int, double, const VectorXd&) -> bool {
        ++call_count; return false;
    };
    MinimizeOptions opts;
    opts.quiet = true;
    opts.rhoend = 1e-4;
    opts.callback = cb;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "no-terminate callback success");
    t.check_close(result.x(0), 1.0, 1e-3, "x[0] close to 1");
    t.check_close(result.x(1), 2.5, 1e-3, "x[1] close to 2.5");
    t.check(call_count > 1, "callback called multiple times");
}

void test_rhoend_without_rhobeg(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    MinimizeOptions opts;
    opts.quiet = true;
    opts.rhoend = 4e-4;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "rhoend only success");
    t.check_close(result.x(0), 1.0, 1e-3, "x[0] close to 1");
    t.check_close(result.x(1), 2.5, 1e-3, "x[1] close to 2.5");
}

void test_rhobeg_without_rhoend(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    MinimizeOptions opts;
    opts.quiet = true;
    opts.rhobeg = -1;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "rhobeg only success");
    t.check_close(result.x(0), 1.0, 1e-3, "x[0] close to 1");
    t.check_close(result.x(1), 2.5, 1e-3, "x[1] close to 2.5");
}

void test_eta2_without_eta1(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    MinimizeOptions opts;
    opts.quiet = true;
    opts.eta2 = 0.7;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "eta2 only success");
    t.check_close(result.x(0), 1.0, 1e-3, "x[0] close to 1");
    t.check_close(result.x(1), 2.5, 1e-3, "x[1] close to 2.5");
}

void test_eta2_out_of_range(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    MinimizeOptions opts;
    opts.quiet = true;
    opts.eta2 = 1.7;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "eta2 out of range success");
    t.check_close(result.x(0), 1.0, 1e-3, "x[0] close to 1");
    t.check_close(result.x(1), 2.5, 1e-3, "x[1] close to 2.5");
}

void test_eta1_out_of_range(TestRunner& t) {
    VectorXd x0(2); x0 << 5.0, 5.0;
    MinimizeOptions opts;
    opts.quiet = true;
    opts.eta1 = 1.7;
    auto result = minimize(obj, x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "eta1 out of range success");
    t.check_close(result.x(0), 1.0, 1e-3, "x[0] close to 1");
    t.check_close(result.x(1), 2.5, 1e-3, "x[1] close to 2.5");
}

void test_method_validation(TestRunner& t) {
    VectorXd x0(2); x0 << 0.0, 0.0;

    bool threw = false;
    try { minimize(obj, x0, "invalid"); }
    catch (const std::invalid_argument&) { threw = true; }
    t.check(threw, "invalid method throws");

    threw = false;
    try { minimize(obj, x0, "newuoa"); }
    catch (const std::invalid_argument&) { threw = true; }
    t.check(threw, "newuoa throws (only cobyla implemented)");
}

void test_minimize_constraint_violation(TestRunner& t) {
    auto combined_nlc = [](const VectorXd& x) -> VectorXd {
        VectorXd r(2);
        r(0) = x(0) - 4;
        r(1) = 5 - x(0);
        return r;
    };
    VectorXd lb2(2);
    lb2 << -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity();
    VectorXd ub2(2); ub2 << 0, 0;
    NonlinearConstraint nlc(combined_nlc, lb2, ub2);
    auto nlc_func = transform_constraint_function(nlc);

    VectorXd x0(1); x0 << 0;
    MinimizeOptions opts;
    opts.quiet = true;
    auto result = minimize([](const VectorXd& x) { return x(0); },
                           x0, "cobyla", nullptr, nullptr, &nlc_func, opts);
    t.check(result.maxcv > 0.1, "constraint violation > 0.1");
}

void test_scalar(TestRunner& t) {
    VectorXd x0(1); x0 << 5;
    MinimizeOptions opts;
    opts.quiet = true;
    auto result = minimize([](const VectorXd& x) { return x(0) * x(0); },
                           x0, "cobyla", nullptr, nullptr, nullptr, opts);
    t.check(result.success, "scalar success");
    t.check_close(result.x(0), 0.0, 1e-3, "x close to 0");
    t.check_close(result.fun, 0.0, 1e-3, "f close to 0");
}

int main() {
    std::cout << "=== test_miscellaneous: Options and Edge Cases ===" << std::endl;
    TestRunner t;

    test_callback_terminate(t);
    test_callback_no_terminate(t);
    test_rhoend_without_rhobeg(t);
    test_rhobeg_without_rhoend(t);
    test_eta2_without_eta1(t);
    test_eta2_out_of_range(t);
    test_eta1_out_of_range(t);
    test_method_validation(t);
    test_minimize_constraint_violation(t);
    test_scalar(t);

    return t.done();
}
