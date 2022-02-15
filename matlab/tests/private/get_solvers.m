function spaths = get_solvers(solvers, compile_flag)

% All possible solvers.
known_solvers = {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'};

% Parse the input.
if nargin == 0
    solvers = known_solvers;
    compile_flag = true;
elseif nargin == 1
    compile_flag = 'true';
    if isa(solvers, 'logical') && numel(solvers) == 1
        compile_flag = solvers;
        solvers = known_solvers;
    end
end
if isa(solvers, 'char') || isa(solvers, 'string')
    solvers = {solvers};
end

% Directories.
test_dir = fileparts(fileparts(mfilename('fullpath'))); % Parent of the directory where this .m file resides.
neupdfo_dir = fileparts(fileparts(test_dir));
callstack = dbstack;
invoker = callstack(2).name;  % The function that calls this function.
if strcmp(invoker, 'verify')
    pdfo_dir = fullfile(neupdfo_dir, 'OPDFO');
else
    pdfo_dir = fullfile(neupdfo_dir, 'PDFO');
end

% Get the solver paths, and add them to the MATLAB path.
solver_dir = fullfile(pdfo_dir, 'matlab', 'interfaces');
solvern_dir = fullfile(neupdfo_dir, 'matlab', 'interfaces');
spaths={solver_dir, solvern_dir};
for ip = 1 : length(spaths)
    addpath(spaths{ip});
end

% Compile the solvers if needed.
if compile_flag
    for is = 1 : length(solvers)
        if strcmp(solvers{is}(end), 'n')
            solver = solvers{is}(1:end-1);
        else
            solver = solvers{is};
        end

        mexopt = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mexopt.classical = false;  % Should be removed when the classical variant is ready
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mexopt.debug = strcmp(invoker, 'verify');

        cd(neupdfo_dir);
        setup([solver, 'n'], mexopt);

        cd(pdfo_dir);
        setup(solver, mexopt)

        cd(test_dir);
    end
end
