function success = set_compiler_options(compiler_options)
%Configure the compiler options by editing the mexopts files. It is hacky!!!

if ~isunix || ismac
    error('Configuration of compiler options supports only Linux.')
end

success = false;

mfilepath = fileparts(mfilename('fullpath')); % The directory containing this setup script

% Modify mexopts files in `config_dir` after making backups.
config_dir = fullfile(matlabroot, 'bin', 'glnxa64', 'mexopts');
config_files = {dir(fullfile(config_dir, 'gfortran*.xml')).name};

try
    fileattrib(config_dir, '+w');
catch
end
if ~iswritable(config_dir)
    warning('The directory %s is not writable. Compile options not set.', config_dir);
    return;
end

for ifile = 1 : length(config_files)
    cfile = fullfile(config_dir, config_files{ifile});

    try
        fileattrib(cfile, '+w')
    catch
    end
    if ~iswritable(cfile)
        warning('The file %s is not writable. Compile options not set.', cfile);
        return;
    end

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

% If MEX has been set up before, then the configuration is already written in the following file.
% We delete it so that MEX will be reconfigured according to `config_files`.
% Why not making a backup for it? Because the existing version may be generated using the modified
% `config_files` when this script was called the last time (N.B.: `mex_setup_file` does not exist
% in a fresh installation of MATLAB), in which case it would be wrong to store this copy and
% restore `mex_setup_file` using it.
mex_setup_file = fullfile(prefdir, ['mex_FORTRAN_', computer('arch'), '.xml']);
if exist(mex_setup_file, 'file')

    try
        fileattrib(prefdir, '+w');
    catch
    end
    if ~iswritable(prefdir)
        warning('The directory %s is not writable. The modified compile options may not take effect.', prefdir);
    end

    try
        fileattrib(mex_setup_file, '+w');
    catch
    end
    if ~iswritable(mex_setup_file)
        warning('The file %s is not writable. The modified compile options may not take effect.', mex_setup_file);
    end

    delete(mex_setup_file);
end

fprintf('\nCompiler options set to \n\n%s\n\n', compiler_options);
success = true;

return
