#ifndef PRIMA_CPP_PREPROC_HPP
#define PRIMA_CPP_PREPROC_HPP

/*
This is a module that preprocesses the inputs.

Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
*/

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "consts.hpp"
#include "present.hpp"

namespace prima {

struct PreprocResult {
    int iprint;
    int maxfun;
    int maxhist;
    double ftarget;
    double rhobeg;
    double rhoend;
    int npt;
    int maxfilt;
    double ctol;
    double cweight;
    double eta1;
    double eta2;
    double gamma1;
    double gamma2;
    Eigen::VectorXd x0;
};

/*
This subroutine preprocesses the inputs. It does nothing to the inputs that are valid.
*/
inline PreprocResult preproc(
    const std::string& solver, int num_vars, int iprint, int maxfun, int maxhist,
    double ftarget, double rhobeg, double rhoend,
    std::optional<int> num_constraints = std::nullopt,
    std::optional<int> npt = std::nullopt,
    std::optional<int> maxfilt = std::nullopt,
    std::optional<double> ctol = std::nullopt,
    std::optional<double> cweight = std::nullopt,
    std::optional<double> eta1 = std::nullopt,
    std::optional<double> eta2 = std::nullopt,
    std::optional<double> gamma1 = std::nullopt,
    std::optional<double> gamma2 = std::nullopt,
    std::optional<bool> is_constrained = std::nullopt,
    std::optional<bool> honour_x0 = std::nullopt,
    std::optional<Eigen::VectorXd> xl = std::nullopt,
    std::optional<Eigen::VectorXd> xu = std::nullopt,
    std::optional<Eigen::VectorXd> x0 = std::nullopt
) {
    // Read num_constraints, if necessary
    int nlc = (num_constraints.has_value() && solver == "cobyla") ? num_constraints.value() : 0;

    // Decide whether the problem is truly constrained
    bool is_con = is_constrained.value_or(nlc > 0);

    // Validate IPRINT
    if (std::abs(iprint) > 3) iprint = IPRINT_DEFAULT;

    // Validate MAXFUN
    double min_maxfun;
    if (solver == "uobyqa") min_maxfun = (num_vars + 1) * (num_vars + 2) / 2 + 1;
    else if (solver == "cobyla") min_maxfun = num_vars + 2;
    else min_maxfun = num_vars + 3;
    if (maxfun < min_maxfun) maxfun = static_cast<int>(min_maxfun);

    // Validate MAXHIST: MAXHIST > MAXFUN is never needed
    if (maxhist < 0) maxhist = maxfun;
    maxhist = std::min(maxhist, maxfun);

    // Validate FTARGET
    if (std::isnan(ftarget)) ftarget = FTARGET_DEFAULT;

    // Validate NPT
    auto npt_val = npt;
    if (npt_val.has_value() && (solver == "newuoa" || solver == "bobyqa" || solver == "lincoa")) {
        if (npt_val.value() < num_vars + 2 || npt_val.value() > std::min(maxfun - 1, ((num_vars + 2) * (num_vars + 1)) / 2)) {
            npt_val = std::min(maxfun - 1, 2 * num_vars + 1);
        }
    }

    // Validate MAXFILT
    auto maxfilt_val = maxfilt;
    if (maxfilt_val.has_value() && (solver == "lincoa" || solver == "cobyla")) {
        int mf_in = maxfilt_val.value();
        if (maxfilt_val.value() < 1) maxfilt_val = MAXFILT_DEFAULT;
        else maxfilt_val = std::max(MIN_MAXFILT, maxfilt_val.value());
        // Further revise MAXFILT according to MAXHISTMEM
        double unit_memo;
        if (solver == "lincoa") unit_memo = (num_vars + 2) * sizeof(double);
        else if (solver == "cobyla") unit_memo = (nlc + num_vars + 2) * sizeof(double);
        else unit_memo = 1;
        // We cannot simply set MAXFILT = MIN(MAXFILT, MAXHISTMEM/...) without type considerations
        if (maxfilt_val.value() > MAXHISTMEM / unit_memo)
            maxfilt_val = static_cast<int>(MAXHISTMEM / unit_memo);
        maxfilt_val = std::min(maxfun, std::max(MIN_MAXFILT, maxfilt_val.value()));
    }

    // Validate ETA1 and ETA2
    double eta1_local = eta1.value_or(ETA1_DEFAULT);
    double eta2_local = eta2.value_or(ETA2_DEFAULT);

    // When the difference between ETA1 and ETA2 is tiny, force them to equal
    if (eta1.has_value() && eta2.has_value()) {
        if (std::abs(eta1.value() - eta2.value()) < 100 * EPS * std::max(std::abs(eta1.value()), 1.0))
            eta2_local = eta1.value();
    }

    if (eta1.has_value()) {
        if (std::isnan(eta1.value())) eta1_local = ETA1_DEFAULT;
        else if (eta1.value() < 0 || eta1.value() >= 1) {
            // Take ETA2 into account if it has a valid value
            if (eta2.has_value() && eta2.value() > 0 && eta2.value() <= 1)
                eta1_local = std::max(EPS, eta2.value() / 7.0);
            else eta1_local = ETA1_DEFAULT;
        } else eta1_local = eta1.value();
    }

    if (eta2.has_value()) {
        if (std::isnan(eta2.value())) eta2_local = ETA2_DEFAULT;
        else if (eta2.value() < eta1_local || eta2.value() > 1)
            eta2_local = (eta1_local + 2) / 3.0;
        else eta2_local = eta2.value();
    }

    // Validate GAMMA1 and GAMMA2
    double gamma1_local = gamma1.value_or(GAMMA1_DEFAULT);
    double gamma2_local = gamma2.value_or(GAMMA2_DEFAULT);

    if (gamma1.has_value()) {
        if (std::isnan(gamma1.value())) gamma1_local = GAMMA1_DEFAULT;
        else if (gamma1.value() <= 0 || gamma1.value() >= 1) gamma1_local = GAMMA1_DEFAULT;
        else gamma1_local = gamma1.value();
    }

    if (gamma2.has_value()) {
        if (std::isnan(gamma2.value())) gamma2_local = GAMMA2_DEFAULT;
        else if (gamma2.value() < 1 || std::isinf(gamma2.value())) gamma2_local = GAMMA2_DEFAULT;
        else gamma2_local = gamma2.value();
    }

    // Validate RHOBEG and RHOEND

    // When the data is passed from interfaces (e.g., MEX) to Fortran code, RHOBEG and RHOEND
    // may change a bit. If we set RHOEND = RHOBEG in the interfaces, it may happen that
    // RHOEND > RHOBEG, which is invalid. Force them to equal when the difference is tiny.
    if (std::abs(rhobeg - rhoend) < 100 * EPS * std::max(std::abs(rhobeg), 1.0))
        rhoend = rhobeg;

    // Revise the default values for RHOBEG/RHOEND according to the solver
    double rhobeg_default = (solver == "bobyqa")
        ? std::max(EPS, std::min(RHOBEG_DEFAULT, (xu.value() - xl.value()).minCoeff() / 4.0))
        : RHOBEG_DEFAULT;
    double rhoend_default = (solver == "bobyqa")
        ? std::max(EPS, std::min(0.1 * rhobeg_default, RHOEND_DEFAULT))
        : RHOEND_DEFAULT;

    if (solver == "bobyqa") {
        // Do NOT merge this IF into the one below! Otherwise, XU and XL may be
        // accessed even if the solver is not BOBYQA.
        // Do NOT make this revision if RHOBEG is not positive or not finite,
        // because otherwise RHOBEG will get a huge value when XU or XL contains
        // huge values that indicate unbounded variables.
        if (rhobeg > (xu.value() - xl.value()).minCoeff() / 2) {
            rhobeg = (xu.value() - xl.value()).minCoeff() / 4.0;
        }
    }
    if (rhobeg <= 0 || std::isnan(rhobeg) || std::isinf(rhobeg)) {
        // Take RHOEND into account if it has a valid value. We do not do
        // this if the solver is BOBYQA, which requires that RHOBEG <= (XU-XL)/2.
        if (std::isfinite(rhoend) && rhoend > 0 && solver != "bobyqa")
            rhobeg = std::max(10 * rhoend, rhobeg_default);
        else rhobeg = rhobeg_default;
    }

    if (rhoend <= 0 || rhobeg < rhoend || std::isnan(rhoend) || std::isinf(rhoend))
        rhoend = std::max(EPS, std::min(0.1 * rhobeg, rhoend_default));

    // For BOBYQA, revise X0 or RHOBEG so that the distance between X0 and
    // the inactive bounds is at least RHOBEG. If HONOUR_X0 == true, revise
    // RHOBEG if needed; otherwise, revise X0 if needed.
    Eigen::VectorXd x0_local = x0.value_or(Eigen::VectorXd());
    if (honour_x0.has_value() && solver == "bobyqa") {
        if (honour_x0.value()) {
            for (int i = 0; i < num_vars; ++i) {
                if (std::isfinite(xl.value()[i]) && std::abs(x0_local[i] - xl.value()[i]) <= EPS * std::max(1.0, std::abs(xl.value()[i])))
                    x0_local[i] = xl.value()[i];
                if (std::isfinite(xu.value()[i]) && std::abs(x0_local[i] - xu.value()[i]) <= EPS * std::max(1.0, std::abs(xu.value()[i])))
                    x0_local[i] = xu.value()[i];
            }
            double min_dist = rhobeg;
            for (int i = 0; i < num_vars; ++i) {
                if (std::isfinite(xl.value()[i]) && std::abs(x0_local[i] - xl.value()[i]) > EPS)
                    min_dist = std::min(min_dist, x0_local[i] - xl.value()[i]);
                if (std::isfinite(xu.value()[i]) && std::abs(x0_local[i] - xu.value()[i]) > EPS)
                    min_dist = std::min(min_dist, xu.value()[i] - x0_local[i]);
            }
            rhobeg = std::max(EPS, min_dist);
        }
    }

    // Validate CTOL (it can be 0)
    if (ctol.has_value()) {
        double c = ctol.value();
        if (std::isnan(c) || c < 0) ctol = CTOL_DEFAULT;
        else ctol = c;
    }

    // Validate CWEIGHT (it can be +Inf)
    if (cweight.has_value()) {
        double cw = cweight.value();
        if (std::isnan(cw) || cw < 0) cweight = CWEIGHT_DEFAULT;
        else cweight = cw;
    }

    return {
        iprint, maxfun, maxhist, ftarget, rhobeg, rhoend,
        npt_val.value_or(0),
        maxfilt_val.value_or(0),
        ctol.value_or(CTOL_DEFAULT),
        cweight.value_or(CWEIGHT_DEFAULT),
        eta1_local, eta2_local, gamma1_local, gamma2_local,
        x0_local
    };
}

} // namespace prima

#endif // PRIMA_CPP_PREPROC_HPP
