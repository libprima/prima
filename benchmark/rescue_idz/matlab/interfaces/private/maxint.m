function mi = maxint()
%MAXINT returns the largest dimension supported. It should be smaller than the largest value of
% mwSize and mwIndex in MEX functions. Otherwise, segmentation faults may occur. For MATLAB R2022a
% on a 64-bit Linux system, this largest value is 2^63 - 1 > 10^18. We set maxint to 0.99 times the
% largest int32 integer, which is about 10^9 and should be safe.
% N.B.: After all, the algorithms are not designed for large problems.

mi = floor(0.99*intmax('int32'));

return
