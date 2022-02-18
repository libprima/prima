function mexname = get_mexname(solver, precision, debug_flag, variant)
%GET_MEXNAME returns the name of the mexified `solver` according to `precision`, `debug_flag`, and
% `variant`.

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

cpwd = fileparts(mfilename('fullpath'));  % The directory where this file resides.

% `precision_list` is a cell array containing all the precisions of the Fortran solvers.
% N.B.: we may call this function during the setup of the package, when `all_precisions.m` does not exist.
allprec_file = fullfile(cpwd, 'all_precisions.m');
if exist(allprec_file, 'file')
    precision_list = all_precisions();
else  % `all_precisions.m` has not been created yet.
    precision_list = all_precisions_possible();
end

% `variant_list` is a cell array containing all the variants of the Fortran solvers.
% N.B.: we may call this function during the setup of the package, when `all_variants.m` does not exist.
allvar_file = fullfile(cpwd, 'all_variants.m');
if exist(allvar_file, 'file')
    variant_list = all_variants();
else  % `all_variants.m` has not been created yet.
    variant_list = all_variants_possible();
end

if nargin ~= 2 && nargin ~= 4
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid number of inputs.', funname);
elseif ~ischarstr(solver) || ~((nargin == 2 && strcmp(solver, 'gethuge')) || (nargin == 4 && ismember(solver, all_solvers())))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid solver received', funname);
elseif ~(ischarstr(precision) && ismember(precision, precision_list) && ~isempty(precision))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid precision received', funname);
elseif nargin == 4 && ~islogicalscalar(debug_flag)
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid debugging flag received', funname);
elseif nargin == 4 && ~(ischarstr(variant) && ismember(variant, variant_list) && ~isempty(precision))
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
