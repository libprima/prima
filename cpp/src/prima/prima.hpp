#ifndef PRIMA_CPP_PRIMA_HPP
#define PRIMA_CPP_PRIMA_HPP

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "infos.hpp"
#include "bounds.hpp"
#include "project.hpp"
#include "linalg.hpp"
#include "linear_constraints.hpp"
#include "nonlinear_constraints.hpp"
#include "evaluate.hpp"
#include "message.hpp"
#include "cobyla/cobyla.hpp"

namespace prima {

struct MinimizeOptions {
    bool quiet = true;
    double rhobeg = std::numeric_limits<double>::quiet_NaN();
    double rhoend = std::numeric_limits<double>::quiet_NaN();
    int maxfun = 0;
    int iprint = 0;
    double ftarget = -std::numeric_limits<double>::infinity();
    double ctol = std::numeric_limits<double>::quiet_NaN();
    double cweight = CWEIGHT_DEFAULT;
    double eta1 = std::numeric_limits<double>::quiet_NaN();
    double eta2 = std::numeric_limits<double>::quiet_NaN();
    double gamma1 = std::numeric_limits<double>::quiet_NaN();
    double gamma2 = std::numeric_limits<double>::quiet_NaN();
    int maxfilt = 2000;
    int maxhist = 0;
    std::function<bool(const Eigen::VectorXd&, double, int, int, double,
                        const Eigen::VectorXd&)> callback = nullptr;
};

struct MinimizeResult {
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

namespace detail {

inline LinearConstraint adjust_linear_constraint_for_fixed(
    const LinearConstraint& lc, const Eigen::VectorXi& fixed_mask,
    const Eigen::VectorXd& fixed_vals, int n_total
) {
    int n_free = n_total - fixed_mask.sum();
    int m = lc.lb.size();
    Eigen::VectorXd contrib = Eigen::VectorXd::Zero(m);
    int f_idx = 0;
    for (int i = 0; i < n_total; ++i)
        if (fixed_mask(i)) { contrib += lc.A.col(i) * fixed_vals(f_idx); ++f_idx; }
    Eigen::MatrixXd new_A(m, n_free);
    int col = 0;
    for (int i = 0; i < n_total; ++i)
        if (!fixed_mask(i)) { new_A.col(col) = lc.A.col(i); ++col; }
    return LinearConstraint(new_A, lc.lb - contrib, lc.ub - contrib);
}

inline NonlinearConstraintFunction wrap_nlc_for_fixed(
    const NonlinearConstraintFunction& nlc, const Eigen::VectorXi& fm,
    const Eigen::VectorXd& fv, int n_total
) {
    return [nlc, fm, fv, n_total](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd newx(n_total);
        int fi = 0, fri = 0;
        for (int i = 0; i < n_total; ++i)
            newx(i) = fm(i) ? fv(fi++) : x(fri++);
        return nlc(newx);
    };
}

inline Eigen::VectorXi compute_fixed_mask(
    const Eigen::VectorXd& lb, const Eigen::VectorXd& ub, double tol
) {
    int n = lb.size();
    Eigen::VectorXi m(n);
    for (int i = 0; i < n; ++i)
        m(i) = (lb(i) <= ub(i) && ub(i) <= lb(i) + tol) ? 1 : 0;
    return m;
}

template <typename T>
std::optional<T> ptr_to_opt(const T* p) {
    if (p) return *p;
    return std::nullopt;
}

inline double compute_maxcv(
    const Eigen::VectorXd& x,
    const Eigen::MatrixXd* A_ineq, const Eigen::VectorXd* b_ineq,
    const Eigen::MatrixXd* A_eq, const Eigen::VectorXd* b_eq,
    const Eigen::VectorXd* nlconstr
) {
    double m = 0;
    if (A_ineq && A_ineq->size() > 0) { double v = (*A_ineq * x - *b_ineq).maxCoeff(); if (v > m) m = v; }
    if (A_eq && A_eq->size() > 0) { double v = (*A_eq * x - *b_eq).cwiseAbs().maxCoeff(); if (v > m) m = v; }
    if (nlconstr && nlconstr->size() > 0) { double v = nlconstr->maxCoeff(); if (v > m) m = v; }
    return m;
}

inline void extract_free_variables(
    const Eigen::VectorXd& full, const Eigen::VectorXi& mask, Eigen::VectorXd& free
) {
    int j = 0;
    for (int i = 0; i < full.size(); ++i)
        if (!mask(i)) free(j++) = full(i);
}

} // namespace detail

inline MinimizeResult minimize(
    std::function<double(const Eigen::VectorXd&)> fun,
    const Eigen::VectorXd& x0,
    const std::string& method = "",
    const Bounds* bounds = nullptr,
    const LinearConstraint* linear_constraint = nullptr,
    const NonlinearConstraintFunction* nonlinear_constraint_function = nullptr,
    MinimizeOptions options = MinimizeOptions{}
) {
    using namespace Eigen;
    using namespace detail;

    MinimizeResult result;
    int lenx0 = x0.size();
    std::string algo = method;
    bool quiet = options.quiet;

    if (algo.empty()) algo = "cobyla";
    else {
        for (auto& c : algo) c = std::tolower(c);
        if (algo != "cobyla")
            throw std::invalid_argument("Only COBYLA is implemented");
    }

    auto [lb, ub] = process_bounds(bounds, lenx0);
    double arr_tol = get_arrays_tol(lb, ub);
    VectorXi fixed_mask = compute_fixed_mask(lb, ub, arr_tol);
    int n_fixed = fixed_mask.sum();
    bool any_fixed = n_fixed > 0;
    bool all_fixed = (n_fixed == lenx0);

    std::function<double(const VectorXd&)> active_fun(std::move(fun));

    if (all_fixed) {
        VectorXd fv(n_fixed);
        for (int i = 0; i < lenx0; ++i)
            fv(i) = std::max(lb(i), std::min(ub(i), 0.5 * (lb(i) + ub(i))));
        double fval = active_fun(fv);
        VectorXd nlconstr;
        if (nonlinear_constraint_function) nlconstr = (*nonlinear_constraint_function)(fv);
        MatrixXd A_eq, A_ineq; VectorXd b_eq, b_ineq;
        if (linear_constraint) {
            auto sep = separate_LC_into_eq_and_ineq(*linear_constraint);
            if (sep.A_eq.size() > 0) { A_eq = sep.A_eq; b_eq = sep.b_eq; }
            if (sep.A_ineq.size() > 0) { A_ineq = sep.A_ineq; b_ineq = sep.b_ineq; }
        }
        double maxcv = compute_maxcv(fv,
            A_ineq.size() > 0 ? &A_ineq : nullptr, b_ineq.size() > 0 ? &b_ineq : nullptr,
            A_eq.size() > 0 ? &A_eq : nullptr, b_eq.size() > 0 ? &b_eq : nullptr,
            nlconstr.size() > 0 ? &nlconstr : nullptr);
        result.x = fv; result.success = true; result.status = FIXED_SUCCESS;
        result.message = "All variables were fixed by the provided bounds.";
        result.fun = fval; result.nfev = 1; result.maxcv = maxcv;
        result.nlconstr = nlconstr; result.method = "cobyla";
        return result;
    }

    VectorXd x0_adj = x0, lb_adj = lb, ub_adj = ub, fixed_vals;
    const LinearConstraint* active_lc = linear_constraint;
    LinearConstraint adjusted_lc;
    NonlinearConstraintFunction active_nlc_fun;
    const NonlinearConstraintFunction* active_nlc = nonlinear_constraint_function;

    if (any_fixed) {
        fixed_vals = VectorXd(n_fixed);
        {
            int j = 0;
            for (int i = 0; i < lenx0; ++i)
                if (fixed_mask(i))
                    fixed_vals(j++) = std::max(lb(i), std::min(ub(i), 0.5 * (lb(i) + ub(i))));
        }
        int n_free = lenx0 - n_fixed;
        x0_adj = VectorXd(n_free); lb_adj = VectorXd(n_free); ub_adj = VectorXd(n_free);
        extract_free_variables(x0, fixed_mask, x0_adj);
        extract_free_variables(lb, fixed_mask, lb_adj);
        extract_free_variables(ub, fixed_mask, ub_adj);

        VectorXi fm = fixed_mask; VectorXd fv = fixed_vals;
        int ntot = lenx0;
        auto orig = std::move(active_fun);
        active_fun = [orig, fm, fv, ntot](const VectorXd& x) -> double {
            VectorXd newx(ntot); int fi = 0, fri = 0;
            for (int i = 0; i < ntot; ++i)
                newx(i) = fm(i) ? fv(fi++) : x(fri++);
            return orig(newx);
        };
        if (linear_constraint) {
            adjusted_lc = adjust_linear_constraint_for_fixed(*linear_constraint, fixed_mask, fixed_vals, lenx0);
            active_lc = &adjusted_lc;
        }
        if (nonlinear_constraint_function) {
            active_nlc_fun = wrap_nlc_for_fixed(*nonlinear_constraint_function, fixed_mask, fixed_vals, lenx0);
            active_nlc = &active_nlc_fun;
        }
    }

    if (active_nlc == nullptr) {
        auto proj = project(x0_adj, lb_adj, ub_adj, active_lc);
        x0_adj = proj.x;
    }

    // Separate linear constraints
    MatrixXd A_eq, A_ineq; VectorXd b_eq, b_ineq;
    if (active_lc) {
        auto sep = separate_LC_into_eq_and_ineq(*active_lc);
        if (sep.A_eq.rows() > 0) { A_eq = sep.A_eq; b_eq = sep.b_eq; }
        if (sep.A_ineq.rows() > 0) { A_ineq = sep.A_ineq; b_ineq = sep.b_ineq; }
    }

    // Evaluate f0 and nlconstr0
    double f0_val = active_fun(x0_adj);
    VectorXd nlconstr0;
    int m_nlcon = 0;
    if (active_nlc) {
        nlconstr0 = (*active_nlc)(x0_adj);
        m_nlcon = nlconstr0.size();
    }

    // Build calcfc like pyprima does: returns (f, nlconstr)
    // Capture active_fun and active_nlc by value to avoid dangling references
    auto fun_copy = active_fun;
    auto nlc_copy = active_nlc ? std::make_shared<NonlinearConstraintFunction>(*active_nlc) : nullptr;
    Calcfc calcfc_inner;
    if (nlc_copy) {
        calcfc_inner = [fun_copy, nlc_copy](const VectorXd& x) -> std::pair<double, VectorXd> {
            return {fun_copy(x), (*nlc_copy)(x)};
        };
    } else {
        calcfc_inner = [fun_copy](const VectorXd& x) -> std::pair<double, VectorXd> {
            return {fun_copy(x), VectorXd(0)};
        };
    }

    // Wrap bounds + linear constraints into AMAT/BVEC via get_lincon
    auto [amat_opt, bvec_opt] = get_lincon(
        ptr_to_opt(A_eq.size() > 0 ? &A_eq : nullptr),
        ptr_to_opt(A_ineq.size() > 0 ? &A_ineq : nullptr),
        ptr_to_opt(b_eq.size() > 0 ? &b_eq : nullptr),
        ptr_to_opt(b_ineq.size() > 0 ? &b_ineq : nullptr),
        ptr_to_opt(&lb_adj),
        ptr_to_opt(&ub_adj)
    );

    auto cobyla_result = cobyla(
        calcfc_inner, m_nlcon, x0_adj,
        amat_opt, bvec_opt, std::nullopt, std::nullopt,
        std::nullopt, std::nullopt,
        f0_val, nlconstr0.size() > 0 ? std::optional<VectorXd>(nlconstr0) : std::nullopt,
        std::isnan(options.rhobeg) ? std::nullopt : std::optional<double>(options.rhobeg),
        std::isnan(options.rhoend) ? std::nullopt : std::optional<double>(options.rhoend),
        options.ftarget,
        std::isnan(options.ctol) ? CTOL_DEFAULT : options.ctol,
        options.cweight,
        options.maxfun > 0 ? std::optional<int>(options.maxfun) : std::nullopt,
        options.iprint,
        std::isnan(options.eta1) ? std::nullopt : std::optional<double>(options.eta1),
        std::isnan(options.eta2) ? std::nullopt : std::optional<double>(options.eta2),
        GAMMA1_DEFAULT, GAMMA2_DEFAULT,
        options.maxhist > 0 ? std::optional<int>(options.maxhist) : std::nullopt,
        options.maxfilt,
        options.callback
    );

    result.x = cobyla_result.x;
    result.fun = cobyla_result.f;
    result.nlconstr = cobyla_result.constr;
    result.maxcv = cobyla_result.cstrv;
    result.nfev = cobyla_result.nf;
    result.status = cobyla_result.info;
    result.success = (cobyla_result.info == SMALL_TR_RADIUS || cobyla_result.info == FTARGET_ACHIEVED);
    result.message = get_info_string("COBYLA", cobyla_result.info);
    result.method = "cobyla";

    if (any_fixed) {
        VectorXd full_x(lenx0);
        int fi = 0, fri = 0;
        for (int i = 0; i < lenx0; ++i)
            full_x(i) = fixed_mask(i) ? fixed_vals(fi++) : result.x(fri++);
        result.x = full_x;
    }

    return result;
}

} // namespace prima

#endif // PRIMA_CPP_PRIMA_HPP
