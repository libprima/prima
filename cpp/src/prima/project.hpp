#ifndef PRIMA_CPP_COMMON_HPP
#define PRIMA_CPP_COMMON_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include "bounds.hpp"
#include "linear_constraints.hpp"
#include "infos.hpp"

namespace prima {

struct OptimizeResult {
    Eigen::VectorXd x;
    bool success = false;
    int status = 0;
    std::string message;
    double fun = 0.0;
    int nfev = 0;
    double maxcv = 0.0;
    Eigen::VectorXd nlconstr;
    std::string method;
};

inline OptimizeResult project(
    const Eigen::VectorXd& x0,
    const Eigen::VectorXd& lb,
    const Eigen::VectorXd& ub,
    const LinearConstraint* linear_constraint
) {
    const double max_con = 1e20;
    const double eps = std::numeric_limits<double>::epsilon();

    OptimizeResult result;

    if (linear_constraint == nullptr) {
        Eigen::VectorXd x_proj = x0.cwiseMax(lb).cwiseMin(ub);
        result.x = x_proj;
        return result;
    }

    const Eigen::MatrixXd& A = linear_constraint->A;
    const Eigen::VectorXd& lc_lb = linear_constraint->lb;
    const Eigen::VectorXd& lc_ub = linear_constraint->ub;

    bool all_eq = true;
    for (Eigen::Index i = 0; i < lc_lb.size(); ++i) {
        if (std::abs(lc_ub(i) - lc_lb(i)) > eps) {
            all_eq = false;
            break;
        }
    }
    bool bounds_trivial = (lb.maxCoeff() <= -max_con) && (ub.minCoeff() >= max_con);

    if (all_eq && bounds_trivial) {
        Eigen::VectorXd b = 0.5 * (lc_lb + lc_ub);
        Eigen::VectorXd rhs = b - A * x0;
#if EIGEN_VERSION_AT_LEAST(5,0,0)
        Eigen::VectorXd xi = A.bdcSvd<Eigen::ComputeThinU | Eigen::ComputeThinV>().solve(rhs);
#else
        Eigen::VectorXd xi = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
#endif
        Eigen::VectorXd x_proj = (x0 + xi).cwiseMax(lb).cwiseMin(ub);
        result.x = x_proj;
        return result;
    }

    if (linear_constraint != nullptr) {
        Eigen::VectorXd ax = A * x0;
        bool feasible = true;

        for (Eigen::Index i = 0; i < lc_lb.size(); ++i) {
            if (ax(i) < lc_lb(i) - eps || ax(i) > lc_ub(i) + eps) {
                feasible = false;
                break;
            }
        }
        for (Eigen::Index i = 0; i < x0.size(); ++i) {
            if (x0(i) < lb(i) - eps || x0(i) > ub(i) + eps) {
                feasible = false;
                break;
            }
        }

        if (feasible) {
            result.x = x0;
            return result;
        }

        Eigen::VectorXd x_proj = x0.cwiseMax(lb).cwiseMin(ub);

        Eigen::MatrixXd Aeq;
        Eigen::VectorXd beq;
        bool has_eq = false;
        {
            std::vector<Eigen::Index> eq_idx;
            for (Eigen::Index i = 0; i < lc_lb.size(); ++i) {
                if (std::abs(lc_ub(i) - lc_lb(i)) <= eps) {
                    eq_idx.push_back(i);
                }
            }
            if (!eq_idx.empty()) {
                Aeq = Eigen::MatrixXd(eq_idx.size(), A.cols());
                beq = Eigen::VectorXd(eq_idx.size());
                for (std::size_t j = 0; j < eq_idx.size(); ++j) {
                    Aeq.row(j) = A.row(eq_idx[j]);
                    beq(j) = 0.5 * (lc_lb(eq_idx[j]) + lc_ub(eq_idx[j]));
                }
                has_eq = true;
            }
        }

        if (has_eq) {
            Eigen::VectorXd rhs = beq - Aeq * x_proj;
#if EIGEN_VERSION_AT_LEAST(5,0,0)
            Eigen::VectorXd step = Aeq.bdcSvd<Eigen::ComputeThinU | Eigen::ComputeThinV>().solve(rhs);
#else
            Eigen::VectorXd step = Aeq.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
#endif
            x_proj = (x_proj + step).cwiseMax(lb).cwiseMin(ub);
        }

        result.x = x_proj;
        return result;
    }

    result.x = x0;
    return result;
}

} // namespace prima

#endif // PRIMA_CPP_COMMON_HPP
