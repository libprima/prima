function mex_solvers(solver)

% Directories.
test_dir = fileparts(fileparts(mfilename('fullpath'))); % Directory where this .m file resides.
neupdfo_dir = fileparts(fileparts(test_dir));

callstack = dbstack;
invoker = callstack(2).name;  % The function that calls this function.
if strcmp(invoker, 'verify')
    pdfo_dir = fullfile(neupdfo_dir, 'OPDFO');
else
    pdfo_dir = fullfile(neupdfo_dir, 'PDFO');
end

% Compile the solvers.
mexopt = struct();
mexopt.debug = strcmp(invoker, 'verify');
cd(neupdfo_dir);
setup([solver, 'n'], mexopt);
cd(pdfo_dir);
setup(solver, mexopt)
cd(test_dir);
