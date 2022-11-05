function [solver_list, options, action, wrong_input] = parse_input(argin)
%PARSE_INPUT parses the input to the setup script.

% Compilation options.
options = [];
action_list = {'compile', 'uninstall', 'clean', 'path'};
action = 'compile';
wrong_input = false;

% Start the parsing to set `input_string` and `options`.
input_string = 'ALL';  % Default value for `input_string`.
if length(argin) > 2
    fprintf('\nSetup accepts at most two inputs.\n\n');
    wrong_input = true;
elseif length(argin) == 1
    if ischarstr(argin{1})
        input_string = argin{1};
    elseif isa(argin{1}, 'struct') || isempty(argin{1})
        options = argin{1};
    else
        fprintf('\nThe input to setup should be a string and/or a structure.\n\n');
        wrong_input = true;
    end
elseif length(argin) == 2
    if (ischarstr(argin{1})) && (isa(argin{2}, 'struct') || isempty(argin{2}))
        input_string = argin{1};
        options = argin{2};
    elseif (ischarstr(argin{2})) && (isa(argin{1}, 'struct') || isempty(argin{1}))
        input_string = argin{2};
        options = argin{1};
    else
        fprintf('\nThe input to setup should be a string and/or a structure.\n\n');
        wrong_input = true;
    end
end

% Cast input_string to a character array in case it is a MATLAB string.
input_string = lower(char(input_string));

% Set `options` to struct() if it is empty.
if isempty(options)
    options = struct();
end

% Parse `input_string` to set `action` and `solver_list`.
if ismember(input_string, action_list)
    solver_list = {};
    action = input_string;
    if ~isempty(options) && ~isempty(fieldnames(options))
        fprintf('\nOptions are ignored since an action (''%s'') is specified.\n', input_string)
    end
else
    solver = input_string;
    % Decide which solver(s) to compile.
    if ismember(solver, all_solvers())
        solver_list = {solver};
    elseif strcmpi(solver, 'ALL')
        solver_list = all_solvers();
    else
        fprintf('Unknown solver ''%s'' to compile.\n\n', solver);
        solver_list = {};
        wrong_input = true;
    end
end

assert(isempty(solver_list) || strcmp(action, 'compile'), '`solver_list` should be empty unless `action` == ''compile''');
assert(isempty(setdiff(solver_list, all_solvers())), '`solver_list` is a subset of `all_solvers()`');
assert(ismember(action, action_list), '`action` belongs to `action_list`');

return
