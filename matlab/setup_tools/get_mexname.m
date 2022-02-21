function mexname = get_mexname(solver, precision, debug_flag, variant, mexdir)
%GET_MEXNAME returns the name of the mexified `solver` according to `precision`, `debug_flag`,
% `variant`, and `mexdir`.
% N.B.:
% 1. `get_mexname` accepts 2, 4, or 5 arguments:
%    get_mexname('gethuge', precision)
%    get_mexname(solver, precision, debug_flag, variant),
%    get_mexname(solver, precision, debug_flag, variant, mexdir),
%    where `solver` is a member of `all_solvers()` in the last two cases.
% 2. `get_mexname` can be called during setup or runtime. During setup, `get_mexname` decides the
%    name of the MEX file to compile; during runtime, it decides the name of MEX file to call.
% 3. When `solver` is 'gethuge', the returns of `get_mexname` are the same during setup and runtime.
% 4. When `solver` is not 'gethuge', `get_mexname` has 4 inputs if it is called during setup and
%    5 inputs if it is called during runtime.
% 3. In general, when `solver` is not 'gethuge', `mexname` will contain the character returned by
%    `dbgstr(debug_flag)`. However, when variant == classical, `mexname` will contain `dbgstr(false)`
%    regardless of `debug_flag`; during runtime , `mexname` will contain either `dbgstr(debug_flag)`
%    or `dbgstr(false)`, depending on the availability of the corresponding MEX file under `mexdir`.

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

solver_list = [all_solvers(), 'gethuge'];

if nargin ~= 2 && nargin ~= 4 && nargin ~= 5
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid number of inputs.', funname);
elseif ~ischarstr(solver) || ~((nargin == 2 && strcmp(solver, 'gethuge')) || ((nargin == 4 || nargin == 5) && ismember(solver, solver_list)))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid solver received', funname);
end

if nargin == 5 && (~ischarstr(mexdir) || isempty(mexdir))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid mexdir received', funname);
end

% `precision_list` is a cell array containing all the precisions of the Fortran solvers.
% `variant_list` is a cell array containing all the variants of the Fortran solvers.
if nargin == 4
    precision_list = all_precisions_possible();
    variant_list = all_variants_possible();
else
    precision_list = all_precisions();
    variant_list = all_variants();
end

if ~(ischarstr(precision) && ismember(precision, precision_list))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid precision received', funname);
elseif (nargin == 4 || nargin == 5) && ~islogicalscalar(debug_flag)
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid debugging flag received', funname);
elseif (nargin == 4 || nargin == 5) && ~(ischarstr(variant) && ismember(variant, variant_list))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid variant received', funname);
end

% Start the real business
if strcmp(solver, 'gethuge')
    mexname = [solver, '_', precision(1)];  % `mexname` is independent of other inputs, if any.
else
    if strcmp(variant, 'classical')
        % The support for the classical variant is limited. We do not provide a debugging version.
        debug_flag = false;
    end
    mexname = [solver, '_', precision(1), dbgstr(debug_flag), variant(1)];
end

if nargin == 5 && ~exist(fullfile(mexdir, mexname), 'file')
    mexname = [solver, '_', precision(1), dbgstr(false), variant(1)];
end

return
