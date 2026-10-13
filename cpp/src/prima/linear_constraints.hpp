#ifndef PRIMA_CPP_LINEAR_CONSTRAINTS_HPP
#define PRIMA_CPP_LINEAR_CONSTRAINTS_HPP

/*
This module handles linear constraints of the form lb <= A*x <= ub.

Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
*/

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <vector>

namespace prima {

// Represents a set of linear constraints: lb <= A*x <= ub
struct LinearConstraint {
    Eigen::MatrixXd A;
    Eigen::VectorXd lb;
    Eigen::VectorXd ub;

    LinearConstraint() = default;

    LinearConstraint(const Eigen::MatrixXd& A_, const Eigen::VectorXd& lb_, const Eigen::VectorXd& ub_)
        : A(A_), lb(lb_), ub(ub_) {}
};

// Combine multiple LinearConstraint objects by vertically stacking
// their A matrices and corresponding lb/ub vectors
inline LinearConstraint combine_multiple_linear_constraints(
    const std::vector<LinearConstraint>& constraints
) {
    if (constraints.empty()) {
        return LinearConstraint();
    }

    Eigen::MatrixXd full_A = constraints[0].A;
    Eigen::VectorXd full_lb = constraints[0].lb;
    Eigen::VectorXd full_ub = constraints[0].ub;

    for (std::size_t i = 1; i < constraints.size(); ++i) {
        Eigen::MatrixXd new_A(full_A.rows() + constraints[i].A.rows(), full_A.cols());
        new_A << full_A, constraints[i].A;
        full_A = new_A;

        Eigen::VectorXd new_lb(full_lb.size() + constraints[i].lb.size());
        new_lb << full_lb, constraints[i].lb;
        full_lb = new_lb;

        Eigen::VectorXd new_ub(full_ub.size() + constraints[i].ub.size());
        new_ub << full_ub, constraints[i].ub;
        full_ub = new_ub;
    }

    return LinearConstraint(full_A, full_lb, full_ub);
}

// Result of separating linear constraints into equality and inequality parts
struct SeparatedLinearConstraints {
    Eigen::MatrixXd A_eq;
    Eigen::VectorXd b_eq;
    Eigen::MatrixXd A_ineq;
    Eigen::VectorXd b_ineq;
};

/*
Separate a LinearConstraint (lb <= A*x <= ub) into:
- Equality constraints: A_eq * x = b_eq (where lb == ub within tolerance)
- Inequality constraints: A_ineq * x <= b_ineq (where lb > -inf or ub < inf)

Lower bounds are converted to the form -A[i,:] * x <= -lb[i].
Upper bounds are converted to the form A[i,:] * x <= ub[i].
*/
inline SeparatedLinearConstraints separate_LC_into_eq_and_ineq(
    const LinearConstraint& linear_constraint
) {
    const double eps = std::numeric_limits<double>::epsilon();
    const double tol = 2.0 * eps;

    const Eigen::MatrixXd& A = linear_constraint.A;
    const Eigen::VectorXd& lb = linear_constraint.lb;
    const Eigen::VectorXd& ub = linear_constraint.ub;

    std::vector<Eigen::Index> eq_indices;
    std::vector<Eigen::Index> ineq_lb_indices;
    std::vector<Eigen::Index> ineq_ub_indices;

    for (Eigen::Index i = 0; i < lb.size(); ++i) {
        if (ub(i) <= lb(i) + tol) {
            eq_indices.push_back(i);
        } else {
            if (lb(i) > -std::numeric_limits<double>::infinity()) {
                ineq_lb_indices.push_back(i);
            }
            if (ub(i) < std::numeric_limits<double>::infinity()) {
                ineq_ub_indices.push_back(i);
            }
        }
    }

    SeparatedLinearConstraints result;

    if (!eq_indices.empty()) {
        result.A_eq = Eigen::MatrixXd(eq_indices.size(), A.cols());
        result.b_eq = Eigen::VectorXd(eq_indices.size());
        for (std::size_t j = 0; j < eq_indices.size(); ++j) {
            result.A_eq.row(j) = A.row(eq_indices[j]);
            result.b_eq(j) = 0.5 * (lb(eq_indices[j]) + ub(eq_indices[j]));
        }
    }

    Eigen::MatrixXd A_ineq_lb, A_ineq_ub;
    Eigen::VectorXd b_ineq_lb, b_ineq_ub;
    bool has_lb = false, has_ub = false;

    if (!ineq_lb_indices.empty()) {
        A_ineq_lb = Eigen::MatrixXd(ineq_lb_indices.size(), A.cols());
        b_ineq_lb = Eigen::VectorXd(ineq_lb_indices.size());
        for (std::size_t j = 0; j < ineq_lb_indices.size(); ++j) {
            A_ineq_lb.row(j) = -A.row(ineq_lb_indices[j]);
            b_ineq_lb(j) = -lb(ineq_lb_indices[j]);
        }
        has_lb = true;
    }

    if (!ineq_ub_indices.empty()) {
        A_ineq_ub = Eigen::MatrixXd(ineq_ub_indices.size(), A.cols());
        b_ineq_ub = Eigen::VectorXd(ineq_ub_indices.size());
        for (std::size_t j = 0; j < ineq_ub_indices.size(); ++j) {
            A_ineq_ub.row(j) = A.row(ineq_ub_indices[j]);
            b_ineq_ub(j) = ub(ineq_ub_indices[j]);
        }
        has_ub = true;
    }

    int total_ineq = (has_lb ? A_ineq_lb.rows() : 0) + (has_ub ? A_ineq_ub.rows() : 0);
    if (total_ineq > 0) {
        result.A_ineq = Eigen::MatrixXd(total_ineq, A.cols());
        result.b_ineq = Eigen::VectorXd(total_ineq);
        int row = 0;
        if (has_lb) {
            result.A_ineq.block(0, 0, A_ineq_lb.rows(), A.cols()) = A_ineq_lb;
            result.b_ineq.segment(0, A_ineq_lb.rows()) = b_ineq_lb;
            row += A_ineq_lb.rows();
        }
        if (has_ub) {
            result.A_ineq.block(row, 0, A_ineq_ub.rows(), A.cols()) = A_ineq_ub;
            result.b_ineq.segment(row, A_ineq_ub.rows()) = b_ineq_ub;
        }
    }

    return result;
}

} // namespace prima

#endif // PRIMA_CPP_LINEAR_CONSTRAINTS_HPP
