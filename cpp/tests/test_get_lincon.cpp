#include <Eigen/Core>
#include <iostream>

#include "prima/bounds.hpp"
#include "prima/linear_constraints.hpp"
#include "prima/cobyla/cobyla.hpp"
#include "test_common.hpp"

using namespace prima;
using namespace Eigen;
using prima::test::TestRunner;

void test_get_lincon(TestRunner& t) {
    MatrixXd Aeq(2, 2); Aeq << 1, 2, 3, 4;
    MatrixXd Aineq(2, 2); Aineq << 5, 6, 7, 8;
    VectorXd beq(2); beq << 9, 10;
    VectorXd bineq(2); bineq << 11, 12;
    VectorXd xl(2); xl << 0, -1;
    VectorXd xu(2); xu << 13, 14;

    auto [amat_opt, bvec_opt] = get_lincon(Aeq, Aineq, beq, bineq, xl, xu);
    t.check(amat_opt.has_value(), "amat has value");
    t.check(bvec_opt.has_value(), "bvec has value");

    MatrixXd amat_expected(10, 2);
    amat_expected << -1, 0, 0, -1, 1, 0, 0, 1, -1, -2, -3, -4, 1, 2, 3, 4, 5, 6, 7, 8;
    VectorXd bvec_expected(10);
    bvec_expected << 0, 1, 13, 14, -9, -10, 9, 10, 11, 12;

    t.check(amat_opt->isApprox(amat_expected), "amat matches");
    t.check(bvec_opt->isApprox(bvec_expected), "bvec matches");
}

void test_get_lincon_boundmax(TestRunner& t) {
    MatrixXd Aeq(2, 2); Aeq << 1, 2, 3, 4;
    MatrixXd Aineq(2, 2); Aineq << 5, 6, 7, 8;
    VectorXd beq(2); beq << 9, 10;
    VectorXd bineq(2); bineq << 11, 12;
    VectorXd xl(2); xl << -BOUNDMAX - 1, -1;
    VectorXd xu(2); xu << 13, 14;

    auto [amat_opt, bvec_opt] = get_lincon(Aeq, Aineq, beq, bineq, xl, xu);
    t.check(amat_opt.has_value(), "boundmax amat has value");
    t.check(bvec_opt.has_value(), "boundmax bvec has value");

    MatrixXd amat_expected(9, 2);
    amat_expected << 0, -1, 1, 0, 0, 1, -1, -2, -3, -4, 1, 2, 3, 4, 5, 6, 7, 8;
    VectorXd bvec_expected(9);
    bvec_expected << 1, 13, 14, -9, -10, 9, 10, 11, 12;

    t.check(amat_opt->isApprox(amat_expected), "boundmax amat matches");
    t.check(bvec_opt->isApprox(bvec_expected), "boundmax bvec matches");
}

void test_get_lincon_none(TestRunner& t) {
    auto [amat, bvec] = get_lincon(std::nullopt, std::nullopt,
                                    std::nullopt, std::nullopt,
                                    std::nullopt, std::nullopt);
    t.check(!amat.has_value(), "none amat is nullopt");
    t.check(!bvec.has_value(), "none bvec is nullopt");
}

void test_combine_linear_constraints(TestRunner& t) {
    MatrixXd A1(1, 2); A1 << 1.0, 0.0;
    MatrixXd A2(1, 2); A2 << 0.0, 1.0;
    LinearConstraint lc1(A1, VectorXd::Constant(1, -1.0), VectorXd::Constant(1, 1.0));
    LinearConstraint lc2(A2, VectorXd::Constant(1, -2.0), VectorXd::Constant(1, 2.0));

    auto combined = combine_multiple_linear_constraints({lc1, lc2});
    t.check(combined.A.rows() == 2, "combined A rows == 2");
    t.check(combined.A.cols() == 2, "combined A cols == 2");
    t.check(std::abs(combined.A(0, 0) - 1.0) < 1e-15, "A(0,0) == 1");
    t.check(std::abs(combined.A(1, 1) - 1.0) < 1e-15, "A(1,1) == 1");
    t.check(std::abs(combined.lb(1) + 2.0) < 1e-15, "lb(1) == -2");
    t.check(std::abs(combined.ub(1) - 2.0) < 1e-15, "ub(1) == 2");
}

void test_separate_LC(TestRunner& t) {
    MatrixXd A(2, 2); A << 1.0, 0.0, 0.0, 1.0;
    LinearConstraint lc(A, Vector2d(2.0, -1.0), Vector2d(2.0, 5.0));
    auto sep = separate_LC_into_eq_and_ineq(lc);
    t.check(sep.A_eq.rows() == 1, "eq rows == 1");
    t.check(sep.A_ineq.rows() == 2, "ineq rows == 2");
    t.check(std::abs(sep.b_eq(0) - 2.0) < 1e-15, "b_eq == 2");

    auto sep2 = separate_LC_into_eq_and_ineq(
        LinearConstraint(A, Vector2d(-1.0, -2.0), Vector2d(3.0, 4.0)));
    t.check(sep2.A_eq.rows() == 0, "no eq rows");
    t.check(sep2.A_ineq.rows() == 4, "4 ineq rows");
}

int main() {
    std::cout << "=== test_get_lincon: Linear Constraint Helpers ===" << std::endl;
    TestRunner t;

    test_get_lincon(t);
    test_get_lincon_boundmax(t);
    test_get_lincon_none(t);
    test_combine_linear_constraints(t);
    test_separate_LC(t);

    return t.done();
}
