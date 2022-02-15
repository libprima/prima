function availability = isavailable(directory, precision_or_variant)
%ISAVAILABLE checks whether a precision or variant is available under `directory` for all the solvers
% returned by `all_solvers()`.
% N.B.: We assume that 'double' and `modern` are always available.

solver_list = all_solvers();

% N.B.: Do NOT call `all_precisions()` to decide `precision_list`! Because `isavailable` is called
% during the setup of the package, when `all_precisions.m` does not exist. Indeed, `isavailable` is
% invoked by `create_all_precisions`, which will create `all_precisions.m`. The same for `variant_list`.
precision_list = {'double', 'single', 'quadruple'};
variant_list = {'modern', 'classical'};


callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if nargin ~= 2 || ~ischarstr(directory) || ~ischarstr(precision_or_variant) || ...
        ~ismember(precision_or_variant, [precision_list, variant_list])
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid input received.', funname);
end

if strcmpi(precision_or_variant, 'double')
    availability = true;
    return;
end

if strcmpi(precision_or_variant, 'modern')
    availability = true;
    return;
end

% `availability` set to true if `directory` contains all solvers of the version specified by
% [precision, debug_flag, variant] with the following `precision`, `debug_flag` and `variant`.

precision = 'double';
if ismember(precision_or_variant, precision_list)  % `precision_or_variant` is a precision.
    precision = precision_or_variant;
end

variant = 'modern';
if ismember(precision_or_variant, variant_list)  % `precision_or_variant` is a variant.
    variant = precision_or_variant;
end

debug_flag = false;  % Default `debug_flag`: non-debugging.

availability = true;
for isol = 1 : length(solver_list)
    mexname = get_mexname(solver_list{isol}, precision, debug_flag, variant);
    if ~exist(fullfile(directory, [mexname, '.', mexext()]), 'file')
        availability = false;
        return
    end
end


% ISAVAILABLE ends
return
