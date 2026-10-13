#ifndef PRIMA_CPP_SELECTX_HPP
#define PRIMA_CPP_SELECTX_HPP

// This module provides subroutines that ensure the returned X is optimal among all the calculated
// points in the sense that no other point achieves both lower function value and lower constraint
// violation at the same time. This module is needed only in the constrained case.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "consts.hpp"
#include "present.hpp"
#include "linalg.hpp"

namespace prima {

// This function compares whether FC1 = (F1, C1) is (strictly) better than FC2 = (F2, C2), which
// basically means that (F1 < F2 and C1 <= C2) or (F1 <= F2 and C1 < C2).
// It takes care of the cases where some of these values are NaN or Inf, even though some cases
// should never happen due to the moderated extreme barrier.
// At return, returns TRUE if and only if (F1, C1) is better than (F2, C2).
// Here, C means constraint violation, which is a nonnegative number.
inline bool isbetter(double f1, double c1, double f2, double c2, double ctol) {
    if (std::isnan(f1) || std::isinf(f1) || std::isnan(c1) || std::isinf(c1)) return false;
    if (std::isnan(f2) || std::isinf(f2) || std::isnan(c2) || std::isinf(c2)) return true;

    // Even though NaN/+Inf should not occur in FC1 or FC2 due to the moderated extreme barrier, for
    // security and robustness, the code below does not make this assumption.
    bool ib = false;
    ib = ib || (f1 < f2 && c1 <= c2);
    ib = ib || (f1 <= f2 && c1 < c2);

    // If C1 <= CTOL and C2 is significantly larger/worse than CTOL, i.e., C2 > MAX(CTOL,CREF),
    // then FC1 is better than FC2 as long as F1 < REALMAX. Normally CREF >= CTOL so MAX(CTOL, CREF)
    // is indeed CREF. However, this may not be true if CTOL > 1E-1*CONSTRMAX.
    double cref = 10 * std::max(EPS, std::min(ctol, 1.0E-2 * CONSTRMAX));  // The MIN avoids overflow.
    ib = ib || (f1 < REALMAX && c1 <= ctol && (c2 > std::max(ctol, cref) || std::isnan(c2)));
    return ib;
}

// This subroutine saves X, F, and CSTRV in XFILT, FFILT, and CFILT (and CONSTR in CONFILT
// if they are present), unless a vector in XFILT[:, :NFILT] is better than X.
// If X is better than some vectors in XFILT[:, :NFILT] then these vectors will be
// removed. If X is not better than any of XFILT[:, :NFILT], but NFILT == MAXFILT,
// then we remove a column from XFILT according to the merit function
// PHI = FFILT + CWEIGHT * max(CFILT - CTOL, 0)
// N.B.:
// 1. Only XFILT[:, :NFILT] and FFILT[:, :NFILT] etc contains valid information,
//    while XFILT[:, NFILT+1:MAXFILT] and FFILT[:, NFILT+1:MAXFILT] etc are not
//    initialized yet.
// 2. We decide whether X is better than another by the ISBETTER function.
inline int savefilt(double cstrv, double ctol, double cweight, double f,
                    const Eigen::VectorXd& x, int nfilt,
                    Eigen::VectorXd& cfilt, Eigen::VectorXd& ffilt,
                    Eigen::MatrixXd& xfilt,
                    const Eigen::VectorXd* constr, Eigen::MatrixXd* confilt) {
    int maxfilt = ffilt.size();

    // Return immediately if any column of XFILT is better than X.
    for (int i = 0; i < nfilt; ++i) {
        if (isbetter(ffilt[i], cfilt[i], f, cstrv, ctol)) return nfilt;
        if (ffilt[i] <= f && cfilt[i] <= cstrv) return nfilt;
    }

    // Decide which columns of XFILT to keep.
    Eigen::VectorXi keep = Eigen::VectorXi::Ones(nfilt);
    for (int i = 0; i < nfilt; ++i) {
        if (isbetter(f, cstrv, ffilt[i], cfilt[i], ctol)) keep[i] = 0;
    }

    // If NFILT == MAXFILT and X is not better than any column of XFILT, then we remove the worst
    // column of XFILT according to the merit function PHI = FFILT + CWEIGHT * MAX(CFILT - CTOL, 0).
    int n_keep = keep.sum();
    if (n_keep == maxfilt) {  // In this case, nfilt == keep.size() == n_keep == maxfilt > 0.
        Eigen::VectorXd cfilt_shifted = (cfilt.array() - ctol).cwiseMax(0);
        Eigen::VectorXd phi;
        if (cweight <= 0) phi = ffilt;
        else if (std::isinf(cweight)) {
            phi = cfilt_shifted;
            // We should not use CFILT here; if MAX(CFILT_SHIFTED) is attained at multiple indices,
            // then we will check FFILT to exhaust the remaining degree of freedom.
        } else {
            phi = ffilt.array().max(-REALMAX).matrix();
            for (int i = 0; i < phi.size(); ++i)
                if (std::isnan(phi[i])) phi[i] = -REALMAX;  // Replace NaN with -REALMAX.
            phi.array() += cweight * cfilt_shifted.array();
        }
        // We select X to maximize PHI. In case there are multiple maximizers, we take the one with
        // the largest CSTRV_SHIFTED; if there are more than one choices, we take the one with the
        // largest F; if there are several candidates, we take the one with the largest CSTRV; if
        // the last comparison still leads to more than one possibilities, then they are equally bad
        // and we choose the first.
        // N.B.:
        // 1. This process is the opposite of selecting KOPT in SELECTX.
        // 2. In finite-precision arithmetic, PHI_1 == PHI_2 and CSTRV_SHIFTED_1 ==
        //    CSTRV_SHIFTED_2 do not ensure that F_1 == F_2!
        double phimax = phi.maxCoeff();
        double cref = 0;
        for (int i = 0; i < nfilt; ++i)
            if (phi[i] >= phimax) cref = std::max(cref, cfilt_shifted[i]);
        double fref = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < nfilt; ++i)
            if (cfilt_shifted[i] >= cref) fref = std::max(fref, ffilt[i]);
        int kworst = 0;
        for (int i = 0; i < nfilt; ++i)
            if (ffilt[i] <= fref && cfilt[i] > cfilt[kworst]) kworst = i;
        if (kworst >= 0 && kworst < nfilt) keep[kworst] = 0;
        n_keep = keep.sum();
    }

    // Keep the good xfilt values and remove all the ones that are strictly worse than the new x.
    int idx = 0;
    for (int i = 0; i < nfilt; ++i) {
        if (keep[i]) {
            xfilt.col(idx) = xfilt.col(i);
            ffilt[idx] = ffilt[i];
            cfilt[idx] = cfilt[i];
            if (confilt && constr) confilt->col(idx) = confilt->col(i);
            ++idx;
        }
    }
    nfilt = idx;

    // Once we have removed all the vectors that are strictly worse than x,
    // we add x to the filter.
    xfilt.col(nfilt) = x;
    ffilt[nfilt] = f;
    cfilt[nfilt] = cstrv;
    if (confilt && constr) confilt->col(nfilt) = *constr;
    ++nfilt;

    return nfilt;
}

// This subroutine selects X according to FHIST and CHIST, which represents (a part of) history
// of F and CSTRV. Normally, FHIST and CHIST are not the full history but only a filter, e.g. ffilt
// and CFILT generated by SAVEFILT. However, we name them as FHIST and CHIST because the [F, CSTRV]
// in a filter should not dominate each other, but this subroutine does NOT assume such a property.
// N.B.: CTOL is the tolerance of the constraint violation (CSTRV). A point is considered feasible if
// its constraint violation is at most CTOL. Note that CTOL is absolute, not relative.
inline int selectx(const Eigen::VectorXd& fhist, const Eigen::VectorXd& chist,
                   double cweight, double ctol) {
    int nhist = fhist.size();

    double fref, cref;
    // We select X among the points with F < FREF and CSTRV < CREF.
    // Do NOT use F <= FREF, because F == FREF (FUNCMAX or REALMAX) may mean F == INF in practice!
    bool has_fc = false;
    for (int i = 0; i < nhist; ++i) {
        if (fhist[i] < FUNCMAX && chist[i] < CONSTRMAX) { has_fc = true; break; }
    }
    if (has_fc) { fref = FUNCMAX; cref = CONSTRMAX; }
    else {
        bool has_fr = false;
        for (int i = 0; i < nhist; ++i)
            if (fhist[i] < REALMAX && chist[i] < CONSTRMAX) { has_fr = true; break; }
        if (has_fr) { fref = REALMAX; cref = CONSTRMAX; }
        else {
            bool has_fc2 = false;
            for (int i = 0; i < nhist; ++i)
                if (fhist[i] < FUNCMAX && chist[i] < REALMAX) { has_fc2 = true; break; }
            if (has_fc2) { fref = FUNCMAX; cref = REALMAX; }
            else { fref = REALMAX; cref = REALMAX; }
        }
    }

    int kopt = nhist - 1;
    bool any_valid = false;
    for (int i = 0; i < nhist; ++i)
        if (fhist[i] < fref && chist[i] < cref) { any_valid = true; break; }

    if (any_valid) {
        // Shift the constraint violations by ctol, so that cstrv <= ctol is regarded as no violation.
        Eigen::VectorXd chist_shifted = (chist.array() - ctol).cwiseMax(0);
        // cmin is the minimal shift constraint violation attained in the history.
        double cmin = std::numeric_limits<double>::infinity();
        for (int i = 0; i < nhist; ++i)
            if (fhist[i] < fref) cmin = std::min(cmin, chist_shifted[i]);
        // We consider only the points whose shifted constraint violations are at most the cref below.
        // N.B.: Without taking std::max(EPS, .), cref would be 0 if cmin = 0. In that case, asking for
        // cstrv_shift < cref would be WRONG!
        double cref2 = std::max(EPS, 2 * cmin);

        // We use the following phi as our merit function to select X.
        Eigen::VectorXd phi(nhist);
        if (cweight <= 0) phi = fhist;
        else if (std::isinf(cweight)) {
            phi = chist_shifted;
            // We should not use chist here; if MIN(chist_shifted) is attained at multiple indices,
            // then we will check fhist to exhaust the remaining degree of freedom.
        } else {
            phi = fhist.array().max(-REALMAX) + cweight * chist_shifted.array();
            // std::max(fhist, -REALMAX) makes sure that phi will not contain NaN (unless there is a bug).
        }

        // We select X to minimize phi subject to f < fref and cstrv_shift <= cref (see the comments
        // above for the reason of taking "<" and "<=" in these two constraints). In case there are
        // multiple minimizers, we take the one with the least cstrv_shift; if there is more than one
        // choice, we take the one with the least f; if there are several candidates, we take the one
        // with the least cstrv; if the last comparison still leads to more than one possibility, then
        // they are equally good and we choose the first.
        // N.B.:
        // 1. This process is the opposite of selecting kworst in savefilt.
        // 2. In finite-precision arithmetic, phi_1 == phi_2 and cstrv_shift_1 == cstrv_shift_2 do
        //    not ensure that f_1 == f_2!
        double phimin = std::numeric_limits<double>::infinity();
        for (int i = 0; i < nhist; ++i)
            if (fhist[i] < fref && chist_shifted[i] <= cref2)
                phimin = std::min(phimin, phi[i]);

        double cmin2 = std::numeric_limits<double>::infinity();
        for (int i = 0; i < nhist; ++i)
            if (fhist[i] < fref && phi[i] <= phimin)
                cmin2 = std::min(cmin2, chist_shifted[i]);

        double fmin2 = std::numeric_limits<double>::infinity();
        for (int i = 0; i < nhist; ++i)
            if (chist_shifted[i] <= cmin2)
                fmin2 = std::min(fmin2, fhist[i]);

        for (int i = 0; i < nhist; ++i)
            if (fhist[i] <= fmin2 && chist[i] < chist[kopt]) kopt = i;
    }

    return kopt;
}

} // namespace prima

#endif // PRIMA_CPP_SELECTX_HPP
