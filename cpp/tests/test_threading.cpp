#include <Eigen/Core>
#include <cmath>
#include <future>
#include <iostream>
#include <vector>

#include "prima/bounds.hpp"
#include "prima/infos.hpp"
#include "prima/linear_constraints.hpp"
#include "prima/nonlinear_constraints.hpp"
#include "prima/prima.hpp"
#include "test_common.hpp"

using namespace prima;
using namespace Eigen;
using prima::test::TestRunner;

static double rosen(const VectorXd& x) {
    double s = 0.0;
    for (int i = 0; i < x.size() - 1; ++i) {
        s += 100.0 * std::pow(x(i + 1) - x(i) * x(i), 2) + std::pow(1.0 - x(i), 2);
    }
    return s;
}

struct Problem {
    VectorXd x0;
    const Bounds* bounds = nullptr;
    const LinearConstraint* lincon = nullptr;
    const NonlinearConstraintFunction* nlcon = nullptr;
    MinimizeOptions opts;
};

static MinimizeResult run_problem(const Problem& p) {
    return minimize(rosen, p.x0, "cobyla", p.bounds, p.lincon, p.nlcon, p.opts);
}

void test_threading(TestRunner& t) {
    // Problem 1: rosen 10 vars, no constraints
    Problem p1;
    p1.x0 = VectorXd::LinSpaced(10, 3.4, 17.8);
    p1.opts.quiet = true;
    p1.opts.rhoend = 1e-4;
    p1.opts.maxfun = 5000;

    // Problem 2: rosen 7 vars, no constraints
    Problem p2;
    p2.x0 = VectorXd::LinSpaced(7, -178.0, -23.0);
    p2.opts.quiet = true;
    p2.opts.rhoend = 1e-4;
    p2.opts.maxfun = 5000;

    // Problem 3: rosen 7 vars + nonlinear constraint: x[0]^2 - 5 >= 0
    Problem p3;
    p3.x0 = VectorXd::LinSpaced(7, -178.0, -23.0);
    auto nlc_lambda = [](const VectorXd& x) -> VectorXd {
        return VectorXd::Constant(1, x(0) * x(0) - 5.0);
    };
    NonlinearConstraint nlc(nlc_lambda, VectorXd::Constant(1, 0.0),
                            VectorXd::Constant(1, std::numeric_limits<double>::infinity()));
    NonlinearConstraintFunction nlc_func = transform_constraint_function(nlc);
    p3.nlcon = &nlc_func;
    p3.opts.quiet = true;
    p3.opts.rhoend = 1e-4;
    p3.opts.maxfun = 5000;

    // Problem 4: rosen 7 vars + bounds + linear equality: sum(x) = 50
    Problem p4;
    p4.x0 = VectorXd::LinSpaced(7, -178.0, -23.0);
    Bounds bounds(VectorXd::Constant(7, 0.0), VectorXd::Constant(7, 10.0));
    p4.bounds = &bounds;
    MatrixXd A(1, 7); A.setOnes();  // sum constraint
    LinearConstraint lc(A, VectorXd::Constant(1, 50.0), VectorXd::Constant(1, 50.0));
    p4.lincon = &lc;
    p4.opts.quiet = true;
    p4.opts.rhoend = 1e-4;
    p4.opts.maxfun = 5000;

    std::vector<Problem> problems = {p1, p2, p3, p4};

    // Single-threaded run
    std::vector<MinimizeResult> single_results;
    for (const auto& prob : problems) {
        single_results.push_back(run_problem(prob));
    }

    // Multi-threaded run
    std::vector<std::future<MinimizeResult>> futures;
    for (const auto& prob : problems) {
        futures.push_back(std::async(std::launch::async, run_problem, prob));
    }
    std::vector<MinimizeResult> multi_results;
    for (auto& f : futures) {
        multi_results.push_back(f.get());
    }

    // Compare results
    for (int i = 0; i < 4; ++i) {
        std::string tag = "problem" + std::to_string(i + 1);
        t.check((single_results[i].x - multi_results[i].x).norm() < 1e-10,
                tag + " x matches");
        t.check(std::abs(single_results[i].fun - multi_results[i].fun) < 1e-10,
                tag + " f matches");
        t.check(single_results[i].nfev == multi_results[i].nfev,
                tag + " nfev matches");
        if (single_results[i].nlconstr.size() == multi_results[i].nlconstr.size()) {
            t.check((single_results[i].nlconstr - multi_results[i].nlconstr).norm() < 1e-10,
                    tag + " constr matches");
        }
    }
}

int main() {
    std::cout << "=== test_threading: Thread-safety of COBYLA ===" << std::endl;
    TestRunner t;
    test_threading(t);
    return t.done();
}
