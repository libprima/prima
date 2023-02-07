function huge = gethuge(data_type, precision)
%GETHUGE calls returns a huge number according to `data_type` and `precision`.

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

if nargin < 1 || nargin > 2
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid number of inputs.', funname);
else
    if ~ischarstr(data_type)
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid data_type received.', funname);
    end
    if nargin == 2 && ~(ischarstr(precision) && ismember(precision, all_precisions()))
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid precision received.', funname);
    end
end

if nargin < 2
    precision = 'double';
end

% The following values are intended to be the returns of the Fortran intrinsics HUGE, RADIX, and
% MAXEXPONENT corresponding to the given precision. They are the largest finite value, base of the
% model that represents the floating point numbers, and maximal exponent of the model.
% As of 20230207, the following values are consistent with gfortran, ifort, ifx, nagfor, nvfortran,
% Classic flang, AOCC flang, sunf95, and g95.
if strcmpi(precision, 'single')
    hugenum = realmax('single') ;
    radix = 2;
    maxexponent = 128;
else  % Even if precision = 'quadruple'
    hugenum = realmax('double') ;
    radix = 2;
    maxexponent = 1024;
end

% The following values are intended to be consistent with HUGENUM, HUGECON, and HUGEFUN defined in
% the Fortran code.
switch lower(data_type)
case {'real'}
    huge = hugenum;
case {'fun', 'function', 'con', 'constraint'}
    huge = radix^min(100, maxexponent / 2);
otherwise
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid data_type received.', funname);
end

return
