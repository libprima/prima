function mexname = get_mexname(solver, precision, debug_flag, variant)
%GET_MEXNAME returns the name of the mexified `solver` according to `precision`, `debug_flag`, and
% `variant`.

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

if nargin ~= 2 && nargin ~= 4
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid number of inputs.', funname);
elseif ~ischarstr(solver) || ~((nargin == 2 && strcmp(solver, 'gethuge')) || (nargin == 4 && ismember(solver, all_solvers())))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid solver received', funname);
elseif ~(ischarstr(precision) && ismember(precision, all_precisions()) && ~isempty(precision))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid precision received', funname);
elseif nargin == 4 && ~islogicalscalar(debug_flag)
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid debugging flag received', funname);
elseif nargin == 4 && ~(ischarstr(variant) && ismember(variant, all_variants()) && ~isempty(precision))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid variant received', funname);
end

switch solver
case {'gethuge'}
    mexname = [solver, '_', precision(1)];
otherwise
    if (strcmp(variant, 'classical'))
        % The support for the classical variant is limited. We do not provide a debugging version.
        debug_flag = false;
    end
    mexname = [solver, '_', precision(1), dbgstr(debug_flag), variant(1)];
end

return
