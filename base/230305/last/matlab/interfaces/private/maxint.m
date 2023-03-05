function mi = maxint()
%MAXINT returns the largest integer in the mex functions; the factor 0.99 provides a buffer. We do
% not pass any integer larger than maxint to the mexified Fortran code, or errors include SEGFAULT
% may occur. The value of maxint is about 10^9 on a 32-bit platform and 10^18 on a 64-bit one.

mi = floor(0.99*min([gethuge('integer'), gethuge('mwSI')]));

return
