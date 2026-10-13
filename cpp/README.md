## About

This is a C++ translation of [Zaikun Zhang](https://www.zhangzk.net)'s [modern-Fortran reference implementation](https://github.com/libprima/prima/tree/main/fortran)
for Powell's derivative-free optimization solvers, which is available at `fortran/` under the root directory.
It is a faithful translation of the [Python translation](https://github.com/libprima/prima/tree/main/pyprima),
following the same structure, variable names, and algorithm logic to keep maintenance across languages tractable.

Due to [bug-fixes](https://github.com/libprima/prima#bug-fixes) and [improvements](https://github.com/libprima/prima#improvements),
the modern-Fortran reference implementation by [Zaikun Zhang](https://www.zhangzk.net)
behaves differently from the original Fortran 77 implementation by [M. J. D. Powell](https://www.zhangzk.net/powell.html),
even though the algorithms are essentially the same. Therefore, it is important to point out that you are using
PRIMA rather than the original solvers if you want your results to be reproducible.

As of June 2026, only the COBYLA solver is available in this C++ translation.
The other solvers will be translated from the Python/Fortran reference implementations in the future.

## Building

This is a header-only library requiring only Eigen3. To build the tests:

```bash
cmake -S cpp -B build -DEigen3_DIR=/path/to/eigen3/cmake
cmake --build build --target test_minimize_cpp_exe
```

To install:

```bash
cmake --install build --prefix /usr/local
```

After installation, use from another project:

```cmake
find_package(primacpp REQUIRED)
target_link_libraries(myapp PRIVATE prima::primacpp)
```

## Development notes

- Function names, variable names, and file layout follow the Fortran and Python implementations.
  Keep them in sync when making changes.
- Comments are kept minimal compared to the Python/Fortran sources. When the intent is unclear,
  refer to the Python or Fortran reference.
- The library is header-only. All implementation is in `.hpp` files under `src/prima/`.
- The namespace is `prima`; internals go in `prima::detail`.
