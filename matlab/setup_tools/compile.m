function compile(solver_list, mex_options, mexdir, src, classical, common, gateways)

% Clean up the directories before compilation.
% This is important especially if there was previously another compilation with different options.
% Without cleaning-up, the MEX files may be linked with wrong files, which can lead to serious
% errors including Segmentation Fault!
dir_list = {mexdir, src, classical, common, gateways};
for idir = 1 : length(dir_list)
    cellfun(@(filename) delete(filename), list_modo_files(dir_list{idir}));
end

% Name of the file that contains the list of Fortran files. There should be such a file in each
% Fortran source code directory, and the list should indicate the dependence among the files.
filelist = 'ffiles.txt';

% Compilation of the common files. They are shared by all solvers. We compile them only once.
% gateways/debug.F contains debugging subroutines tailored for MEX.
copyfile(fullfile(gateways, 'debug.F'), common);
% ppf.h contains preprocessing directives. Set __DEBUGGING__ according to mex_options.
header_file = fullfile(common, 'ppf.h');
header_file_bak = fullfile(common, 'ppf.h.bak');
copyfile(header_file, header_file_bak);
if ismember('-g', mex_options)
    rep_str(header_file, '#define __DEBUGGING__ 0', '#define __DEBUGGING__ 1');
else
    rep_str(header_file, '#define __DEBUGGING__ 1', '#define __DEBUGGING__ 0');
end
% Common Fortran source files.
common_files = regexp(fileread(fullfile(common, filelist)), '\n', 'split');
common_files = strtrim(common_files(~cellfun(@isempty, common_files)));
common_files = fullfile(common, common_files);
common_files = [common_files, fullfile(gateways, 'fmxapi.F'), fullfile(gateways, 'cbfun.F')];
% The loop below may be written in one line as follows:
%mex(mex_options{:}, '-c', common_files{:}, '-outdir', common);
% But it does not work for some versions of MATLAB. This may be because the compilation above does
% not respect the order of common_files{:}, which is critical due to the dependence among modules.
for icf = 1 : length(common_files)
    mex(mex_options{:}, '-c', common_files{icf}, '-outdir', common);
end
common_obj_files = list_obj_files(common);

% Compilation of function gethuge
gateway = fullfile(gateways, 'gethuge.F');
mexname = 'gethuge';
mex(mex_options{:}, common_obj_files{:}, gateway, '-output', mexname, '-outdir', mexdir);

version_list = {'m', 'c'};  % m - modernized, c - classical
for isol = 1 : length(solver_list)

    solver = solver_list{isol};
    fprintf('Compiling %s ... ', solver);

    gateway = fullfile(gateways, [solver, '_mex.F']);

    for iver = 1 : length(version_list)
        switch version_list{iver}
        case 'm'
            srcdir = fullfile(src, solver);
            mexname = ['f', solver, 'n'];
        case 'c'
            srcdir = fullfile(classical, solver);
            mexname = ['f', solver, 'n_classical'];
        end

        % Clean up the source file directory
        cellfun(@(filename) delete(filename), list_modo_files(srcdir));
        % Compile
        src_files = regexp(fileread(fullfile(srcdir, filelist)), '\n', 'split');
        src_files = strtrim(src_files(~cellfun(@isempty, src_files)));
        src_files = fullfile(srcdir, src_files);
        for isf = 1 : length(src_files)
            mex(mex_options{:}, '-c', src_files{isf}, '-outdir', srcdir);
        end
        obj_files = [common_obj_files, list_obj_files(srcdir)];
        mex(mex_options{:}, obj_files{:}, gateway, '-output', mexname, '-outdir', mexdir);
        % Clean up the source file directory
        cellfun(@(filename) delete(filename), list_modo_files(srcdir));
    end
    fprintf('Done.\n');
end

% Clean up common.
cellfun(@(filename) delete(filename), list_modo_files(common));
% Clean up mexdir.
cellfun(@(filename) delete(filename), list_modo_files(mexdir));

% Restore header_file
if exist(header_file_bak, 'file')
    movefile(header_file_bak, header_file);
end
