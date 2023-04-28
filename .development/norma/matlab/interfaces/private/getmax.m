function maxnum = getmax(data_type, precision)
%GETMAX calls returns a huge number according to `data_type` and `precision`.

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

% The following values are intended to be the returns of the Fortran intrinsics RADIX, HUGE, and
% MAXEXPONENT corresponding to the given precision. They are the largest finite value, base of the
% model that represents the floating point numbers, and maximal exponent of the model.
% As of 20230207, the following values are consistent with gfortran, ifort, ifx, nagfor, nvfortran,
% Classic flang, AOCC flang, sunf95, and g95.
radix = 2;
if strcmpi(precision, 'single')
    maxfloat = realmax('single') ;
    maxexponent = 128;
else  % Even if precision = 'quadruple'
    maxfloat = realmax('double') ;
    maxexponent = 1024;
end

% The following values are intended to be consistent with BOUNDMAX, FUNCMAX, and CONSTRMAX defined
% in the Fortran code.
switch lower(data_type)
case {'real'}
    maxnum = maxfloat;
case {'bound'}
    maxnum = 0.25 * maxfloat;
case {'fun', 'func', 'function', 'con', 'constr', 'constraint'}
    maxnum = radix^min(100, maxexponent / 2);
otherwise
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid data_type received.', funname);
end

return
