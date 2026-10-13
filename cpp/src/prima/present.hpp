#ifndef PRIMA_CPP_PRESENT_HPP
#define PRIMA_CPP_PRESENT_HPP

/*
This is a C++ equivalent of the Fortran 'present' function for optional arguments.
*/

#include <optional>

namespace prima {

// present check for raw pointers (nullptr check)
template <typename T>
inline bool present(const T* x) { return x != nullptr; }

// present check for std::optional
template <typename T>
inline bool present(const std::optional<T>& x) { return x.has_value(); }

} // namespace prima

#endif // PRIMA_CPP_PRESENT_HPP
