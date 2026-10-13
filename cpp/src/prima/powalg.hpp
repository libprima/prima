#ifndef PRIMA_CPP_POWALG_HPP
#define PRIMA_CPP_POWALG_HPP

// This module provides some Powell-style linear algebra procedures.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <tuple>

#include "consts.hpp"
#include "linalg.hpp"

namespace prima {

// This function updates the QR factorization of an MxN matrix A of full column rank, attempting to
// add a new column C to this matrix as the LAST column while maintaining the full-rankness.
// Case 1. If C is not in range(A) (theoretically, it implies N < M), then the new matrix is
//          [A, C] (hstack).
// Case 2. If C is in range(A), then the new matrix is [A(:, :n-1), C] (hstack).
// N.B.:
//   0. Instead of R, this subroutine updates Rdiag, which is diag(R), with a size at most M and at
//      least min(m, n+1). The number is min(m, n+1) rather than min(m, n) as n may be augmented by 1
//      in the function.
//   1. With the two cases specified as above, this function does not need A as an input.
//   2. The function changes only Q[:, nsave:m] (nsave is the original value of n) and
//      R[:, n-1] (n takes the updated value).
//   3. Indeed, when C is in range(A), Powell wrote in comments that "set iOUT to the index of the
//      constraint (here, column of A --- Zaikun) to be deleted, but branch if no suitable index can be
//      found". The idea is to replace a column of A by C so that the new matrix still has full rank
//      (such a column must exist unless C = 0). But his code essentially sets iout=n always. Maybe he
//      found this worked well enough in practice. Meanwhile, Powell's code includes a snippet that can
//      never be reached, which was probably intended to deal with the case that IOUT != n.
inline std::tuple<Eigen::MatrixXd, Eigen::VectorXd, int>
qradd_Rdiag(const Eigen::VectorXd& c, Eigen::MatrixXd& Q,
            const Eigen::VectorXd& Rdiag_in, int n) {
    int m = Q.cols();
    int nsave = n;  // Needed for debugging (only)
    Eigen::VectorXd Rdiag = Rdiag_in;

    // As in Powell's COBYLA, CQ is set to 0 at the positions with CQ being negligible as per ISMINOR.
    // This may not be the best choice if the subroutine is used in other contexts, e.g. LINCOA.
    Eigen::VectorXd cq = matprod(c, Q);
    Eigen::VectorXd ca = c.cwiseAbs();
    Eigen::MatrixXd Qa = Q.cwiseAbs();
    Eigen::VectorXd cqa = matprod(ca, Qa);
    // The line below basically makes an element of cq 0 if adding it to the corresponding element of
    // cqa does not change the latter.
    for (int i = 0; i < cq.size(); ++i) {
        if (isminor(cq[i], cqa[i])) cq[i] = 0;
    }

    // Update Q so that the columns of Q[:, n+1:m] are orthogonal to C. This is done by applying a 2D
    // Givens rotation to Q[:, [k, k+1]] from the right to zero C' @ Q[:, k+1] out for K=n+1, ... m-1.
    // Nothing will be done if n >= m-1.
    for (int k = m - 2; k >= std::max(0, n - 1); --k) {
        // Powell wrote cq[k+1] != 0 instead of abs. The two differ if cq[k+1] is NaN.
        // If we apply the rotation below when cq[k+1] = 0, then cq[k] will get updated to |cq[k]|.
        if (std::abs(cq[k + 1]) > 0) {
            Eigen::MatrixXd G = planerot(Eigen::Vector2d(cq[k], cq[k + 1]));
            Eigen::MatrixXd Qk = Q.block(0, k, Q.rows(), 2) * G.transpose();
            Q.block(0, k, Q.rows(), 2) = Qk;
            cq[k] = hypot(cq[k], cq[k + 1]);
        }
    }

    // Augment n by 1 if C is not in range(A)
    if (n < m) {
        // Powell's condition for the following if: cq[n+1] != 0
        if (std::abs(cq[n]) > EPS * EPS && !isminor(cq[n], cqa[n])) {
            n += 1;
        }
    }

    // Update Rdiag so that Rdiag[n] = cq[n] = dot(c, q[:, n]). Note that N may have been augmented.
    if (n - 1 >= 0 && n - 1 < m) {  // n >= m should not happen unless the input is wrong
        Rdiag[n - 1] = cq[n - 1];
    }

    return {Q, Rdiag, n};
}

// This function updates the QR factorization for an MxN matrix A=Q*R so that the updated Q and
// R form a QR factorization of [A_0, ..., A_{I-1}, A_{I+1}, ..., A_{N-1}, A_I] which is the matrix
// obtained by rearranging columns [I, I+1, ... N-1] of A to [I+1, ..., N-1, I]. Here A is ASSUMED
// TO BE OF FULL COLUMN RANK, Q is a matrix whose columns are orthogonal, and R, which is not
// present, is an upper triangular matrix whose diagonal entries are nonzero. Q and R need not be
// square.
// N.B.:
//   0. Instead of R, this function updates Rdiag, which is diag(R), the size being n.
//   1. With L = Q.cols() = R.rows(), we have M >= L >= N. Most often L = M or N.
//   2. This function changes only Q[:, i:] and Rdiag[i:].
//   3. (NDB 20230919) In Python, i is either icon or nact - 2, whereas in FORTRAN it is either
//      icon or nact - 1.
// Used in COBYLA.
inline std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
qrexc_Rdiag(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Q_in,
            const Eigen::VectorXd& Rdiag_in, int i) {
    // Sizes
    int m = A.rows(), n = A.cols();

    // Preconditions
    // assert(n >= 1 && n <= m);
    // assert(i >= 0 && i < n);
    // assert(Rdiag.size() == n);
    // assert(Q.rows() == m && Q.cols() >= n && Q.cols() <= m);
    // double tol = std::max(1.0E-8, std::min(1.0E-1, 1.0E8 * EPS * m + 1));
    // assert(isorth(Q, tol));  // Costly!

    Eigen::MatrixXd Q = Q_in;
    Eigen::VectorXd Rdiag = Rdiag_in;

    if (i < 0 || i >= n) return {Q, Rdiag};

    // Let R be the upper triangular matrix in the QR factorization, namely R = Q^T * A.
    // For each k, find the Givens rotation G with G * (R[k:k+2, :]) = [hypt, 0], and update
    // Q[:, k:k+2] to Q[:, k:k+2] * G^T. Then R = Q^T * A is an upper triangular matrix as long
    // as A[:, [k, k+1]] is updated to A[:, [k+1, k]]. Indeed, this new upper triangular matrix can
    // be obtained by first updating R[[k, k+1], :] to G * (R[[k, k+1], :]) and then exchanging its
    // columns K and K+1; at the same time, entries k and k+1 of R's diagonal Rdiag become
    // [hypt, -(Rdiag[k+1]/hypt) * Rdiag[k]].
    // After this is done for each k = 0, ..., n-2, we obtain the QR factorization of the matrix
    // that rearranges columns [i, i+1, ... n-1] of A as [i+1, ..., n-1, i].
    // Powell's code, however, is slightly different: before everything, he first exchanged columns
    // k and k+1 of Q (as well as rows k and k+1 of R). This makes sure that the entries of the
    // updated Rdiag are all positive if it is the case for the original Rdiag.
    for (int k = i; k < n - 1; ++k) {
        Eigen::Vector2d gvec(Rdiag[k + 1], inprod(Q.col(k), A.col(k + 1)));
        Eigen::MatrixXd G = planerot(gvec);
        Eigen::MatrixXd cols = Q.middleCols(k + 1, 1);
        Eigen::MatrixXd Qpair(Q.rows(), 2);
        Qpair.col(0) = Q.col(k + 1);
        Qpair.col(1) = Q.col(k);
        Qpair = Qpair * G.transpose();
        Q.col(k) = Qpair.col(0);
        Q.col(k + 1) = Qpair.col(1);
        // Powell's code updates Rdiag in the following way:
        // hypt = sqrt(Rdiag[k+1]^2 + dot(Q[:, k], A[:, k+1])^2)
        // Rdiag[[k, k+1]] = [hypt, (Rdiag[k+1]/hypt) * Rdiag[k]]
        // Note that Rdiag[n-1] inherits all rounding in Rdiag[i:n-1] and Q[:, i:n-1] and hence
        // contains significant errors. Thus we may modify Powell's code to set only
        // Rdiag[k] = hypt here and then calculate Rdiag[n] by an inner product after the loop.
        // Nevertheless, we simply calculate Rdiag from scratch below.
    }

    // Calculate Rdiag(i:n) from scratch
    for (int k = i; k < n - 1; ++k) {
        Rdiag[k] = inprod(Q.col(k), A.col(k + 1));
    }
    Rdiag[n - 1] = inprod(Q.col(n - 1), A.col(i));

    return {Q, Rdiag};
}

} // namespace prima

#endif // PRIMA_CPP_POWALG_HPP
