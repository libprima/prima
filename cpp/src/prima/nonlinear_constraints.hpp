#ifndef PRIMA_CPP_NONLINEAR_CONSTRAINTS_HPP
#define PRIMA_CPP_NONLINEAR_CONSTRAINTS_HPP

/*
This module handles nonlinear constraints of the form lb <= constraint_fun(x) <= ub.

Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
*/

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace prima {

// Represents a nonlinear constraint: lb <= fun(x) <= ub
struct NonlinearConstraint {
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> fun;
    Eigen::VectorXd lb;
    Eigen::VectorXd ub;

    NonlinearConstraint() = default;

    NonlinearConstraint(
        std::function<Eigen::VectorXd(const Eigen::VectorXd&)> fun_,
        const Eigen::VectorXd& lb_,
        const Eigen::VectorXd& ub_
    ) : fun(std::move(fun_)), lb(lb_), ub(ub_) {}
};

using NonlinearConstraintFunction = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

/*
Transform a NonlinearConstraint into a function that returns constraint violations.

For each component i:
  - If lb(i) > -inf: returns lb(i) - values(i) (negative if satisfied)
  - If ub(i) <  inf: returns values(i) - ub(i) (negative if satisfied)

The output is a concatenated vector of all active constraint violations.
If lb or ub has size 1 but values has size > 1, the scalar bound is broadcast.
*/
inline NonlinearConstraintFunction transform_constraint_function(
    const NonlinearConstraint& nlc
) {
    return [nlc](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd values = nlc.fun(x);

        Eigen::VectorXd lb = nlc.lb;
        Eigen::VectorXd ub = nlc.ub;

        // Broadcast scalar bounds to match the output size
        if (lb.size() == 1 && values.size() > 1) {
            lb = Eigen::VectorXd::Constant(values.size(), nlc.lb(0));
        }
        if (ub.size() == 1 && values.size() > 1) {
            ub = Eigen::VectorXd::Constant(values.size(), nlc.ub(0));
        }

        if (values.size() != lb.size()) {
            throw std::invalid_argument(
                "The number of elements in the constraint function's output "
                "does not match the number of elements in the lower bound."
            );
        }
        if (values.size() != ub.size()) {
            throw std::invalid_argument(
                "The number of elements in the constraint function's output "
                "does not match the number of elements in the upper bound."
            );
        }

        // Count result size: one entry per finite bound
        int n_out = 0;
        for (Eigen::Index i = 0; i < lb.size(); ++i) {
            if (lb(i) > -std::numeric_limits<double>::infinity()) ++n_out;
            if (ub(i) < std::numeric_limits<double>::infinity()) ++n_out;
        }

        Eigen::VectorXd result(n_out);
        int idx = 0;
        for (Eigen::Index i = 0; i < lb.size(); ++i) {
            if (lb(i) > -std::numeric_limits<double>::infinity()) {
                result(idx++) = lb(i) - values(i);
            }
        }
        for (Eigen::Index i = 0; i < ub.size(); ++i) {
            if (ub(i) < std::numeric_limits<double>::infinity()) {
                result(idx++) = values(i) - ub(i);
            }
        }

        return result;
    };
}

/*
Combine multiple NonlinearConstraint functions into a single function.
Each constraint is transformed into a violation function, and the results
are concatenated into a single output vector.
*/
inline NonlinearConstraintFunction process_nl_constraints(
    const std::vector<NonlinearConstraint>& nlcs
) {
    std::vector<NonlinearConstraintFunction> functions;
    functions.reserve(nlcs.size());
    for (const auto& nlc : nlcs) {
        functions.push_back(transform_constraint_function(nlc));
    }

    return [functions](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        std::vector<Eigen::VectorXd> all_values;
        int total_size = 0;
        for (const auto& fun : functions) {
            Eigen::VectorXd v = fun(x);
            total_size += v.size();
            all_values.push_back(std::move(v));
        }
        Eigen::VectorXd result(total_size);
        int idx = 0;
        for (const auto& v : all_values) {
            result.segment(idx, v.size()) = v;
            idx += v.size();
        }
        return result;
    };
}

} // namespace prima

#endif // PRIMA_CPP_NONLINEAR_CONSTRAINTS_HPP
