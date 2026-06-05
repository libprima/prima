#ifndef PRIMA_CPP_HISTORY_HPP
#define PRIMA_CPP_HISTORY_HPP

// This module provides subroutines that handle the X/F/C histories of the solver, taking into
// account that MAXHIST may be smaller than NF.
//
// Translated from Zaikun Zhang's modern-Fortran reference implementation in PRIMA.
//
// Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

#include <Eigen/Core>
#include <vector>

namespace prima {

// Save the data values to the history lists.
inline void savehist(int maxhist, const Eigen::VectorXd& x,
                     std::vector<Eigen::VectorXd>& xhist,
                     double f, std::vector<double>& fhist,
                     double cstrv, std::vector<double>& chist,
                     const Eigen::VectorXd& constr, std::vector<Eigen::VectorXd>& conhist) {
    if (static_cast<int>(xhist.size()) < maxhist) {
        xhist.push_back(x);
        fhist.push_back(f);
        chist.push_back(cstrv);
        conhist.push_back(constr);
    } else {
        // This effectively accomplishes what rangehist does in the Fortran implementation.
        xhist.erase(xhist.begin());
        fhist.erase(fhist.begin());
        chist.erase(chist.begin());
        conhist.erase(conhist.begin());
        xhist.push_back(x);
        fhist.push_back(f);
        chist.push_back(cstrv);
        conhist.push_back(constr);
    }
}

} // namespace prima

#endif // PRIMA_CPP_HISTORY_HPP
