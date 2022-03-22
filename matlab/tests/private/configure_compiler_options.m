function configure_compiler_options(compiler_options)
%Configure the compiler options by editing the mexopts files. It is hacky!!!

if ~isunix || ismac
    error('Configuration of compiler options supports only Linux.')
end

mfilepath = fileparts(mfilename('fullpath')); % The directory containing this setup script
setup_tools = fullfile(fileparts(fileparts(mfilepath)), 'setup_tools');
addpath(setup_tools);  % We use `rep_str` from `setup_tools`.

config_dir = fullfile(matlabroot,'bin', 'glnxa64', 'mexopts');
config_files = {dir(fullfile(config_dir, 'gfortran*.xml')).name};
fileattrib(config_dir, '+w');
for ifile = 1 : length(config_files)
    cfile = fullfile(config_dir, config_files{ifile});
    fileattrib(cfile, '+w')

    cfile_orig = fullfile(config_dir, [config_files{ifile}, '.orig']);
    if ~exist(cfile_orig, 'file')
        % This will be true if the script is called for the first time.
        copyfile(cfile, cfile_orig, 'f');
    else
        % Restore the original configuration file before configuring the compiler options.
        copyfile(cfile_orig, cfile, 'f');
    end

    cfile_bak = fullfile(config_dir, [config_files{ifile}, '.bak']);
    if exist(cfile_bak, 'file')
        delete(cfile_bak);
    end
    copyfile(cfile, cfile_bak, 'f');

    % Configure the compiler options by editing `cfile`.
    rep_str(cfile, 'FOPTIMFLAGS="-O2"', ['FOPTIMFLAGS="', compiler_options, '"']);
    rep_str(cfile, 'FDEBUGFLAGS="-g"', ['FDEBUGFLAGS="', compiler_options, '"']);
    rep_str(cfile, 'LDOPTIMFLAGS="-O"', ['LDOPTIMFLAGS="', compiler_options, '"']);
    rep_str(cfile, 'LDDEBUGFLAGS="-g"', ['LDDEBUGFLAGS="', compiler_options, '"']);
end

rmpath(setup_tools);  % Remove `setup_tools` from path since it has finishes its job.

% If MEX has been set up before, then the configuration is already written in the following file.
% We back it up and delete it so that MEX will be reconfigured according to `config_files`.
fileattrib(prefdir, '+w');
mex_setup_file = fullfile(prefdir, ['mex_FORTRAN_', computer('arch'), '.xml']);
mex_setup_file_orig = fullfile(prefdir, ['mex_FORTRAN_', computer('arch'), '.xml.orig']);
mex_setup_file_bak = fullfile(prefdir, ['mex_FORTRAN_', computer('arch'), '.xml.bak']);

if exist(mex_setup_file, 'file') && ~exist(mex_setup_file_orig, 'file')
    % This will be true if the script is called for the first time.
    copyfile(mex_setup_file, mex_setup_file_orig, 'f');
end
if exist(mex_setup_file, 'file')
    movefile(mex_setup_file, mex_setup_file_bak, 'f');
end

return
