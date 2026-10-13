#ifndef PRIMA_CPP_LINALG_HPP
#define PRIMA_CPP_LINALG_HPP

/*
This module provides some basic linear algebra procedures.

Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
*/

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "consts.hpp"

namespace prima {

/*
We use naive implementations of matrix multiplication and other routines for two
reasons:
1. When Fortran is compiled in debug mode, and C++ is using these routines, we
   can get bit for bit identical results as compared to Fortran. This is helpful
   for comparing the two implementations. It will be particularly helpful when porting
   the other implementations like LINCOA, etc.
2. On some problems this algorithm is very sensitive to errors in finite precision
   arithmetic. Switching to naive implementation will slow down the algorithm, but
   may be more stable.
*/
inline bool USE_NAIVE_MATH = false;

/*
This function tests whether x is minor compared to ref. It is used by Powell, e.g., in COBYLA.
In precise arithmetic, isminor(x, ref) is true if and only if x == 0; in floating point
arithmetic, isminor(x, ref) is true if x is 0 or its nonzero value can be attributed to
computer rounding errors according to ref.
Larger sensitivity means the function is more strict/precise, the value 0.1 being due to Powell.

For example:
isminor(1e-20, 1e300) -> true, because in floating point arithmetic 1e-20 cannot be added to
1e300 without being rounded to 1e300.
isminor(1e300, 1e-20) -> false, because in floating point arithmetic adding 1e300 to 1e-20
dominates the latter number.
isminor(3, 4) -> false, because 3 can be added to 4 without being rounded off
*/
inline bool isminor(double x, double ref) {
    double sensitivity = 0.1;
    double refa = std::abs(ref) + sensitivity * std::abs(x);
    double refb = std::abs(ref) + 2 * sensitivity * std::abs(x);
    return std::abs(ref) >= refa || refa >= refb;
}

// Element-wise isminor for vectors
inline Eigen::VectorXi isminor(const Eigen::VectorXd& x, const Eigen::VectorXd& ref) {
    Eigen::VectorXi r(x.size());
    for (Eigen::Index i = 0; i < x.size(); ++i) r[i] = isminor(x[i], ref[i]) ? 1 : 0;
    return r;
}

// Dot product; uses Eigen's dot or naive implementation depending on USE_NAIVE_MATH
inline double inprod(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    if (!USE_NAIVE_MATH) return x.dot(y);
    double r = 0;
    for (Eigen::Index i = 0; i < x.size(); ++i) r += x[i] * y[i];
    return r;
}

// Matrix-vector product for 1x2 shapes (row vector * matrix)
inline Eigen::VectorXd matprod12(const Eigen::VectorXd& x, const Eigen::MatrixXd& y) {
    Eigen::VectorXd r(y.cols());
    for (Eigen::Index i = 0; i < y.cols(); ++i) r[i] = inprod(x, y.col(i));
    return r;
}

// Matrix-vector product for 2x1 shapes (matrix * column vector)
inline Eigen::VectorXd matprod21(const Eigen::MatrixXd& x, const Eigen::VectorXd& y) {
    Eigen::VectorXd r(x.rows());
    r.setZero();
    for (Eigen::Index i = 0; i < x.cols(); ++i) r += x.col(i) * y[i];
    return r;
}

// Matrix-matrix product (naive triple-loop)
inline Eigen::MatrixXd matprod22(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y) {
    Eigen::MatrixXd r(x.rows(), y.cols());
    r.setZero();
    for (Eigen::Index i = 0; i < y.cols(); ++i)
        for (Eigen::Index j = 0; j < x.cols(); ++j)
            r.col(j) += x.col(i) * y(i, j);
    return r;
}

// Dispatches matprod12 or uses Eigen's multiply based on USE_NAIVE_MATH
inline Eigen::VectorXd matprod(const Eigen::VectorXd& x, const Eigen::MatrixXd& y) {
    if (!USE_NAIVE_MATH) return x.transpose() * y;
    return matprod12(x, y);
}

// Dispatches matprod21 or uses Eigen's multiply based on USE_NAIVE_MATH
inline Eigen::VectorXd matprod(const Eigen::MatrixXd& x, const Eigen::VectorXd& y) {
    if (!USE_NAIVE_MATH) return x * y;
    return matprod21(x, y);
}

// Dispatches matprod22 or uses Eigen's multiply based on USE_NAIVE_MATH
inline Eigen::MatrixXd matprod(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y) {
    if (!USE_NAIVE_MATH) return x * y;
    return matprod22(x, y);
}

// Outer product; uses Eigen or naive implementation
inline Eigen::MatrixXd outprod(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    if (!USE_NAIVE_MATH) return x * y.transpose();
    Eigen::MatrixXd r(x.size(), y.size());
    for (Eigen::Index i = 0; i < y.size(); ++i) r.col(i) = x * y[i];
    return r;
}

/*
Solve the linear least squares problem using Powell's method.
If USE_NAIVE_MATH is false, uses Eigen's SVD-based solver.
Otherwise, uses the modified Gram-Schmidt / QR approach as in Powell's implementation,
where columns are processed in reverse order with isminor checks to detect
rank-deficient systems.
*/
inline Eigen::VectorXd lsqr(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                             const Eigen::MatrixXd& Q, const Eigen::VectorXd& Rdiag) {
    int n = A.cols();
    if (n == 0) return Eigen::VectorXd();

    // Use the Powell QR-based solver with precomputed Q and Rdiag (the diagonal of R).
    // This is the same algorithm as the Fortran lsqr_Rdiag, which is the original
    // COBYLA approach. See PRIMA's fortran/common/linalg.f90.
    // In precise arithmetic, the algorithm performs back-substitution on the
    // upper-triangular system R * x = Q^T * b. It forces x[j] = 0 when the
    // computed value can be attributed to computer rounding errors (isminor).
    int rank = std::min(static_cast<int>(A.rows()), n);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd y = b;
    for (int i = rank - 1; i >= 0; --i) {
        double yq = inprod(y, Q.col(i));
        double yqa = inprod(y.cwiseAbs(), Q.col(i).cwiseAbs());
        if (isminor(yq, yqa)) { x[i] = 0; }
        else { x[i] = yq / Rdiag[i]; y -= x[i] * A.col(i); }
    }
    return x;
}

// Compute sqrt(x1^2 + x2^2) with overflow/underflow protection
inline double hypot(double x1, double x2) {
    if (!USE_NAIVE_MATH) return std::hypot(x1, x2);
    if (!std::isfinite(x1)) return std::abs(x1);
    if (!std::isfinite(x2)) return std::abs(x2);
    double a = std::min(std::abs(x1), std::abs(x2));
    double b = std::max(std::abs(x1), std::abs(x2));
    if (a > std::sqrt(REALMIN) && b < std::sqrt(REALMAX / 2.1)) {
        return std::sqrt(a * a + b * b);
    } else if (b > 0) {
        return b * std::sqrt((a / b) * (a / b) + 1);
    }
    return 0;
}

// Euclidean norm; uses Eigen or naive implementation
inline double norm(const Eigen::VectorXd& x) {
    if (!USE_NAIVE_MATH) return x.norm();
    // NOTE: Avoid std::pow! And exponentiation in general!
    double r = 0;
    for (Eigen::Index i = 0; i < x.size(); ++i) r += x[i] * x[i];
    return std::sqrt(r);
}

/*
According to its documentation, Eigen's sum may sometimes do partial pairwise summation.
For our purposes, when comparing, we want don't want to do anything fancy, and we
just want to add things up one at a time.
*/
inline double primasum(const Eigen::VectorXd& x) { return x.sum(); }

inline double primasum(const Eigen::MatrixXd& x) { return x.sum(); }

// Sum along specified axis
inline Eigen::VectorXd primasum(const Eigen::MatrixXd& x, int axis) {
    if (axis == 0) return x.colwise().sum();
    return x.rowwise().sum();
}

// Element-wise square
inline Eigen::MatrixXd primapow2(const Eigen::MatrixXd& x) { return x.cwiseProduct(x); }

/*
Believe it or not, x**2 is not always the same as x*x in some languages.
In Fortran they appear to be identical, but in C++ we always use x*x.
*/
inline Eigen::VectorXd primapow2(const Eigen::VectorXd& x) { return x.cwiseProduct(x); }

/*
As in MATLAB, planerot(x) returns a 2x2 Givens matrix G for x in R2 so that Y = G @ x has Y[1] = 0.
Roughly speaking, G = [[x[0]/R, x[1]/R], [-x[1]/R, x[0]/R]], where R = norm(x).
0. We need to take care of the possibilities of R=0, Inf, NaN, and over/underflow.
1. The G defined above is continuous with respect to X except at 0. Following this definition,
   G = [[sign(x[0]), 0], [0, sign(x[0])]] if x[1] == 0,
   G = [[0, sign(x[1])], [sign(x[1]), 0]] if x[0] == 0
   Yet some implementations ignore the signs, leading to discontinuity and numerical instability.
2. Difference from MATLAB: if x contains NaN or consists of only Inf, MATLAB returns a NaN matrix,
   but we return an identity matrix or a matrix of +/-sqrt(2). We intend to keep G always orthogonal.
*/
inline Eigen::MatrixXd planerot(const Eigen::Vector2d& x) {
    // Define C = X(1) / R and S = X(2) / R with R = HYPOT(X(1), X(2)). Handle Inf/NaN, over/underflow.
    double c, s;
    if (std::isnan(x[0]) || std::isnan(x[1])) {
        // In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
        c = 1; s = 0;
    }
    else if (!std::isfinite(x[0]) && !std::isfinite(x[1])) {
        // In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
        c = 1.0 / std::sqrt(2) * (x[0] >= 0 ? 1 : -1);
        s = 1.0 / std::sqrt(2) * (x[1] >= 0 ? 1 : -1);
    }
    // X(1) == 0 == X(2)
    else if (std::abs(x[0]) <= 0 && std::abs(x[1]) <= 0) { c = 1; s = 0; }
    /*
    N.B.:
    0. With <= instead of <, this case covers X(1) == 0 == X(2), which is treated above separately
       to avoid the confusing SIGN(., 0).
    1. SIGN(A, 0) = ABS(A) in Fortran but sign(0) = 0 in MATLAB, Python, Julia, and R.
    2. Taking SIGN(X(1)) into account ensures the continuity of G with respect to X except at 0.
    */
    else if (std::abs(x[1]) <= EPS * std::abs(x[0])) { c = x[0] >= 0 ? 1 : -1; s = 0; }
    /*
    N.B.: SIGN(A, X) = ABS(A) * sign of X != A * sign of X. Therefore, it is WRONG to define G
    as SIGN(RESHAPE([ZERO, -ONE, ONE, ZERO], [2, 2]), X(2)). This mistake was committed on
    20211206 and took a whole day to debug! NEVER use SIGN on arrays unless you are really sure.
    */
    else if (std::abs(x[0]) <= EPS * std::abs(x[1])) { c = 0; s = x[1] >= 0 ? 1 : -1; }
    else {
        /*
        Here is the normal case. It implements the Givens rotation in a stable & continuous way as in:
        Bindel, D., Demmel, J., Kahan, W., and Marques, O. (2002). On computing Givens rotations
        reliably and efficiently. ACM Transactions on Mathematical Software (TOMS), 28(2), 206-238.
        N.B.: 1. Modern compilers compute SQRT(REALMIN) and SQRT(REALMAX/2.1) at compilation time.
        2. The direct calculation without involving T and U seems to work better; use it if possible.
        */
        if ((std::sqrt(REALMIN) < std::abs(x[0]) && std::abs(x[0]) < std::sqrt(REALMAX / 2.1)) &&
            (std::sqrt(REALMIN) < std::abs(x[1]) && std::abs(x[1]) < std::sqrt(REALMAX / 2.1))) {
            // Do NOT use HYPOTENUSE here; the best implementation for one may be suboptimal for the other
            double r = norm(x);
            c = x[0] / r;
            s = x[1] / r;
        } else if (std::abs(x[0]) > std::abs(x[1])) {
            double t = x[1] / x[0];
            double u = std::max(1.0, std::max(std::abs(t), std::sqrt(1 + t * t))); // MAXVAL: precaution against rounding error.
            u *= (x[0] >= 0 ? 1 : -1);
            c = 1 / u;
            s = t / u;
        } else {
            double t = x[0] / x[1];
            double u = std::max(1.0, std::max(std::abs(t), std::sqrt(1 + t * t))); // MAXVAL: precaution against rounding error.
            u *= (x[1] >= 0 ? 1 : -1);
            c = t / u;
            s = 1 / u;
        }
    }
    Eigen::Matrix2d G;
    G << c, s, -s, c; // G = [c, s; -s, c]
    return G;
}

// Set entries of cvshift to zero where they are minor compared to cvsabs
inline void apply_isminor(Eigen::VectorXd& cvshift, const Eigen::VectorXd& cvsabs) {
    for (Eigen::Index i = 0; i < cvshift.size(); ++i) {
        if (isminor(cvshift[i], cvsabs[i])) cvshift[i] = 0;
    }
}

// Test if A is lower triangular within tolerance
inline bool istril(const Eigen::MatrixXd& A, double tol = 0) {
    return (A.cwiseAbs() - A.triangularView<Eigen::Lower>().toDenseMatrix().cwiseAbs()).lpNorm<Eigen::Infinity>() <= tol;
}

// Test if A is upper triangular within tolerance
inline bool istriu(const Eigen::MatrixXd& A, double tol = 0) {
    return (A.cwiseAbs() - A.triangularView<Eigen::Upper>().toDenseMatrix().cwiseAbs()).lpNorm<Eigen::Infinity>() <= tol;
}

// Compute A^T * B
inline Eigen::MatrixXd matprod_transpose(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
    return A.transpose() * B;
}

/*
Compute the inverse of A. Uses Eigen's inverse by default.
If USE_NAIVE_MATH is true, exploits triangular structure (used in COBYLA)
or falls back to QR-based inversion.
*/
inline Eigen::MatrixXd inv(const Eigen::MatrixXd& A) {
    if (!USE_NAIVE_MATH) return A.inverse();
    int n = A.rows();
    if (istril(A)) {
        // This case is invoked in COBYLA.
        Eigen::MatrixXd R = A.transpose();
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; ++i) {
            B(i, i) = 1.0 / R(i, i);
            if (i > 0) {
                Eigen::MatrixXd Bsub = B.block(0, 0, i, i);
                B.block(0, i, i, 1) = -matprod(Bsub, Eigen::MatrixXd(R.block(0, i, i, 1))) / R(i, i);
            }
        }
        return B.transpose();
    } else if (istriu(A)) {
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; ++i) {
            B(i, i) = 1.0 / A(i, i);
            if (i > 0) {
                Eigen::MatrixXd Bsub = B.block(0, 0, i, i);
                B.block(0, i, i, 1) = -matprod(Bsub, Eigen::MatrixXd(A.block(0, i, i, 1))) / A(i, i);
            }
        }
        return B;
    } else {
        // This is NOT the best algorithm for the inverse, but since the QR subroutine is available ...
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
        Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>().toDenseMatrix();
        Eigen::MatrixXd Q = qr.householderQ();
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, n);
        for (int i = n - 1; i >= 0; --i) {
            Eigen::MatrixXd Bsub = B.rightCols(n - 1 - i);
            Eigen::MatrixXd Rsub = R.block(i + 1, i, n - 1 - i, 1);
            B.col(i) = (Q.col(i) - matprod(Bsub, Rsub)) / R(i, i);
        }
        return B;
    }
}

// QR decomposition with column pivoting
inline std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXi> qr(const Eigen::MatrixXd& A) {
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    Eigen::MatrixXd Q = qr.householderQ();
    Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>().toDenseMatrix();
    Eigen::VectorXi P = qr.colsPermutation().indices();
    return {Q, R, P};
}

// Test whether A = B^{-1} up to tolerance tol
inline bool isinv(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, double tol = 0) {
    int n = A.rows();
    double t = (tol > 0) ? tol : std::min(1e-3, 100.0 * EPS * std::max(n, 1));
    t = std::max({t, t * A.cwiseAbs().maxCoeff(), t * B.cwiseAbs().maxCoeff()});
    return ((matprod(A, B) - Eigen::MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff() <= t) ||
           ((matprod(B, A) - Eigen::MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff() <= t);
}

// Test whether matrix A has orthonormal columns up to tolerance tol
inline bool isorth(const Eigen::MatrixXd& A, double tol = 0) {
    int n = A.cols();
    if (n > A.rows()) return false;
    if (std::isnan(A.cwiseAbs().sum())) return false;
    Eigen::MatrixXd At = A.transpose();
    Eigen::MatrixXd M = matprod(At, A) - Eigen::MatrixXd::Identity(n, n);
    if (tol > 0) return M.cwiseAbs().maxCoeff() <= std::max(tol, tol * A.cwiseAbs().maxCoeff());
    return M.cwiseAbs().maxCoeff() <= 0;
}

/*
Get a relative tolerance for a set of arrays. Borrowed from COBYQA.

Parameters:
  arrays: vector of vectors to get the tolerance for.

Returns:
  Relative tolerance for the set of arrays.

Throws:
  std::invalid_argument if no array is provided.
*/
inline double get_arrays_tol(const std::vector<Eigen::VectorXd>& arrays) {
    if (arrays.empty()) throw std::invalid_argument("At least one array must be provided.");
    Eigen::Index size = 0;
    double weight = 1.0;
    for (const auto& arr : arrays) {
        size = std::max(size, arr.size());
        double m = 0;
        for (Eigen::Index i = 0; i < arr.size(); ++i)
            if (std::isfinite(arr[i])) m = std::max(m, std::abs(arr[i]));
        weight = std::max(weight, m);
    }
    if (size == 0) size = 1;
    return 10.0 * EPS * static_cast<double>(size) * weight;
}

// Convenience overload for two vectors
inline double get_arrays_tol(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    return get_arrays_tol({a, b});
}

} // namespace prima

#endif // PRIMA_CPP_LINALG_HPP
