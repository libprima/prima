function availability = isavailable(directory, precision_or_variant)
%ISAVAILABLE checks whether a precision or variant is available under `directory` for all the solvers
% returned by `all_solvers()`.
% N.B.: We assume that `default_precision` and `default_variant` (see below) are always available.

solver_list = all_solvers();

% N.B.: Do NOT call `all_precisions()` to decide `precision_list`! Because `isavailable` is called
% during the setup of the package (if `setup path` is called), when `all_precisions.m` does not exist.
% Indeed, `isavailable` is invoked by `create_all_precisions`, which will create `all_precisions.m`.
% The same for `variant_list`.
[precision_list, default_precision] = all_precisions_possible();
[variant_list, default_variant] = all_variants_possible();


callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if nargin ~= 2 || ~ischarstr(directory) || ~ischarstr(precision_or_variant) || ...
        ~ismember(precision_or_variant, [precision_list, variant_list])
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid input received.', funname);
end

if strcmpi(precision_or_variant, default_precision)
    availability = true;
    return;
end

if strcmpi(precision_or_variant, default_variant)
    availability = true;
    return;
end

% `availability` set to true if `gethuge` and all the solvers present in `directory` are available
% with the version specified by the following [precision, default_debug_flag, variant].
% N.B.: If no solver is present in `directory`, then `availability` is totally defined by the
% availability of `gethuge`.
precision = default_precision;
if ismember(precision_or_variant, precision_list)  % `precision_or_variant` is a precision.
    precision = precision_or_variant;
end
variant = default_variant;
if ismember(precision_or_variant, variant_list)  % `precision_or_variant` is a variant.
    variant = precision_or_variant;
end
default_debug_flag = false;  % Default debug flag: non-debugging.

% Check the availability of `gethuge` corresponding to `precision` (note that the availability
% of `gethuge` does not depend on `default_debug_flag` or `variant`). If it is unavailable, then set
% `availability` to false and return.
mexname = get_mexname('gethuge', precision);
if ~exist(fullfile(directory, [mexname, '.', mexext()]), 'file')
    availability = false;
    return
end

% If `gethuge` is available corresponding to [precision, default_debug_flag, variant], then we set
% `availability` to true if and only if [precision, default_debug_flag, variant] is available for all the
% solvers that are present under `directory`. A solver is considered present if it is available for
% at least one combination of precision and variant with the default debug flag. The following lines
% decide the presence of the solvers.
solver_present = false(length(solver_list), 1);
for isol = 1 : length(solver_list)
    for iprc = 1 : length(precision_list)
        for ivar = 1 : length(variant_list)
            mexname = get_mexname(solver_list{isol}, precision_list{iprc}, default_debug_flag, variant_list{ivar});
            if exist(fullfile(directory, [mexname, '.', mexext()]), 'file')
                solver_present(isol) = true;
                break
            end
        end
        if solver_present(isol)
            break
        end
    end
end

% Now set `availability` according to all the solvers that are present under `directory`.
availability = any(solver_present);  % Default `availability` to true unless there is no solver present.
if availability
    for isol = 1 : length(solver_list)
        if solver_present(isol)
            mexname = get_mexname(solver_list{isol}, precision, default_debug_flag, variant);
            if ~exist(fullfile(directory, [mexname, '.', mexext()]), 'file')
                availability = false;
                return
            end
        end
    end
end


% ISAVAILABLE ends
return
