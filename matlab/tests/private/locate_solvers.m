function spaths = locate_solvers()

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
solver_dir = fullfile(pdfo_dir, 'matlab', 'interfaces');
addpath(solver_dir);

solvern_dir = fullfile(neupdfo_dir, 'matlab', 'interfaces');
addpath(solvern_dir);

spaths={solver_dir, solvern_dir};
