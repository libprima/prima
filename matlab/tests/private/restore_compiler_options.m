function restore_compiler_options()

if ~isunix || ismac
    error('Configuration of compiler options supports only Linux.')
end

config_dir = fullfile(matlabroot,'bin', 'glnxa64', 'mexopts');
config_files = {dir(fullfile(config_dir, 'gfortran*.xml')).name};
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
                warning('Two different back-ups found:\n\n%s\n%s\n\nRestoration done by\n%s\n\n', cfile_orig, cfile_bak, cfile_orig);
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

mex_setup_file = fullfile(prefdir, ['mex_FORTRAN_', computer('arch'), '.xml']);
mex_setup_file_orig = fullfile(prefdir, ['mex_FORTRAN_', computer('arch'), '.xml.orig']);
mex_setup_file_bak = fullfile(prefdir, ['mex_FORTRAN_', computer('arch'), '.xml.bak']);
if exist(mex_setup_file_orig, 'file')
    copyfile(mex_setup_file_orig, mex_setup_file, 'f');
    if exist(mex_setup_file_bak, 'file')
        difference = [];
        if isunix
            [~, difference] = system(['diff ''', mex_setup_file_orig, ''' ''', mex_setup_file_bak, '''']);
        end
        if ~isempty(difference)
            warning('Two different back-ups found:\n\n%s\n%s\n\nRestoration done by\n%s\n\n', mex_setup_file_orig, mex_setup_file_bak, mex_setup_file_orig);
        else
            delete(mex_setup_file_bak);
        end
    end
elseif exist(mex_setup_file_bak, 'file')
    movefile(mex_setup_file_bak, mex_setup_file, 'f');
else
    delete(mex_setup_file);
    warning('Failed to restore %s.\nIt is deleted and will be regenerated when MEX is set up for the next time.\n', mex_setup_file);
end

return
