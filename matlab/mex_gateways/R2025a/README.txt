 This directory contains temporary MEX gateway to circumvent the MATLAB R2025a bug that it segfaults
 when the Fortran MEX function contains an internal procedure that is passed as an actual argument.

 The MEX files use a module variable to store the function handle, which is essentially a
 global variable and is not thread-safe or recursion-safe.

 See MathWorks Technical Support Case 07931486 and
 https://www.mathworks.com/matlabcentral/answers/2178414-bug-matlab-2025a-segfaults-on-ubuntu-when-handling-fortran-mex-files-with-internal-subroutines
 https://stackoverflow.com/questions/79699706/matlab-2025a-vs-fortran-mex-files-with-internal-subroutines
 https://fortran-lang.discourse.group/t/implementation-of-a-parametrized-objective-function-without-using-module-variables-or-internal-subroutines
 https://stackoverflow.com/questions/79705107/fortran-implementating-a-parametrized-objective-function-without-using-module-v

Zaikun ZHANG (www.zhangzk.net)
July 20, 2025
