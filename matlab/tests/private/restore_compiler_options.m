function success = restore_compiler_options()
%Restore the compiler options by restoring the mexopts files. It is hacky!!!

if ~isunix || ismac
    error('Configuration of compiler options supports only Linux.')
end

success = false;

config_dir = fullfile(matlabroot,'bin', 'glnxa64', 'mexopts');
config_files = {dir(fullfile(config_dir, 'gfortran*.xml')).name};

try
    fileattrib(config_dir, '+w');
catch
end
if ~iswritable(config_dir)
    warning('The directory %s is not writable. Compile options not restored.', config_dir);
    return;
end

for ifile = 1 : length(config_files)
    cfile = fullfile(config_dir, config_files{ifile});
    cfile_orig = fullfile(config_dir, [config_files{ifile}, '.orig']);
    cfile_bak = fullfile(config_dir, [config_files{ifile}, '.bak']);

    if exist(cfile_orig, 'file')
        copyfile(cfile_orig, cfile, 'f');
        if exist(cfile_bak, 'file')
            difference = [];
            if isunix
                [~, difference] = system(['diff ''', cfile_orig, ''' ''', cfile_bak, '''']);
            end
            if ~isempty(difference)
                warning('Two different backups found:\n\n%s\n%s\n\nRestoration done by\n%s\n\n', cfile_orig, cfile_bak, cfile_orig);
            else
                delete(cfile_bak);
            end
        end
    elseif exist(cfile_bak, 'file')
        movefile(cfile_bak, cfile, 'f');
    else
        error('Failed to restore %s', cfile);
    end
end

% Delete `mex_setup_file` so that it will be regenerated next time when `mex -setup` is called.
% Why not restore it using a backup? See the comments in `set_compiler_options` when this file
% was deleted.
mex_setup_file = fullfile(prefdir, ['mex_FORTRAN_', computer('arch'), '.xml']);
if exist(mex_setup_file, 'file')

    try
        fileattrib(prefdir, '+w');
    catch
    end
    if ~iswritable(prefdir)
        warning('The directory %s is not writable. The restored compile options may not take effect.', prefdir);
    end

    try
        fileattrib(mex_setup_file, '+w');
    catch
    end
    if ~iswritable(mex_setup_file)
        warning('The file %s is not writable. The restored compile options may not take effect.', mex_setup_file);
    end

    delete(mex_setup_file);
end

fprintf('\nCompiler options restored successfully.\n\n');
success = true;

return
