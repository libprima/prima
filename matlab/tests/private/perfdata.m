function summary_file = perfdata(solvers, options)

data_dir = options.data_dir;
test_feature = options.test_feature;
feature_and_time = [test_feature, '.', options.time];

% `stamp` is an identifier for the test. N.B.:
% 1. `stamp` does not include `time`, so that it is easy to identify when we reload the data.
% 2. `stamp` does not include `test_feature` (but the file names will), so that we can make
%    `test_feature` appear after `perf_*`, which is convenient when we view the files in terminal
%    (when using tab for auto-completion, we do not need to type `test_feature` if the files in the
%    current directory correspond to only one `test_feature`).
stamp = strcat(strjoin(solvers, '_'), '.', int2str(options.mindim), '_', int2str(options.maxdim), '.', options.type);
options.stamp = strcat(stamp, '.', feature_and_time);

% `outdir` is the directory to contain the figures (.eps, .pdf) and problem list (problem.txt).
outdir = fullfile(data_dir, strcat(stamp, '.', feature_and_time));
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Check whether to reload data or calculate everything from scratch.
matfile = fullfile(data_dir, strcat(stamp, '.perfdata.', test_feature, '.mat'));
if (isfield(options, 'reload') && options.reload == true)
    % Use directly the data of the last computation for the same solvers and dimensions.
    % `stamp_with_wildcard` replaces `options.type` in `stamp` by a wildcard that contains `options.type`.
    % Be careful with the replacement, as `prob.type` may consist of only one letter (e.g., 'u').
    stamp_with_wildcard = regexprep(stamp, ['\.', options.type, '$'], ['\.\*', options.type, '\*']);
    matfile_with_wildcard = strrep(matfile, stamp, stamp_with_wildcard);
    mlist = dir(matfile_with_wildcard);
    if ~isempty(mlist)
        % Take the latest .mat file that complies with the name pattern `matfile_with_wildcard`.
        [~, im] = max([mlist(:).datenum]);
        matfile = fullfile(mlist(im).folder, mlist(im).name);
    else
        error('Fail to load data, because no .mat file found with the name pattern \n\n%s', matfile_with_wildcard);
    end
    load(matfile, 'frec', 'fmin', 'pdim', 'plist');
    fprintf('\nSucceed in loading data from\n\n%s\n\n', matfile)
    qlist = secup(options);
    fun = @(p) ismember(p, qlist);
    pind = find(cellfun(fun, plist));
    frec = frec(pind, :, :, :);
    fmin = fmin(pind);
    pdim = pdim(pind);
    plist = plist(pind);
else
    if exist(matfile, 'file')
        delete(matfile);
    end
    [frec, fmin, output] = testcu(solvers, options);
    pdim = output.pdim;
    plist = output.plist;
    save(matfile, 'frec', 'fmin', 'pdim', 'plist', '-v7.3');
end

% Record the tested problems in `problems.txt`.
fprob = fullfile(outdir, strcat(stamp, '.', feature_and_time, '.', 'problems.txt'));
fid = fopen(fprob, 'w');
if fid >= 3
    for ip = 1 : length(plist)
        fprintf(fid, "%s\n", plist{ip});
    end
    fclose(fid);
else
    error('\nFail to open file %s.\n', fprob);
end

% Plot the profiles.
ns = length(solvers);
for is = 1:ns
    solvers{is} = regexprep(solvers{is}, '_4test', '');
    solvers{is} = regexprep(solvers{is}, '_classical$', ' (classical)');
    solvers{is} = regexprep(solvers{is}, '_single$', ' (single)');
    solvers{is} = regexprep(solvers{is}, '_quadruple$', ' (quadruple)');
    solvers{is} = regexprep(solvers{is}, '_archiva$', ' (archiva)');
    solvers{is} = regexprep(solvers{is}, '_norma$', ' (norma)');
    %solvers{is} = regexprep(solvers{is}, 'newuoa', 'NEWUOA');
end

prof_options = struct();
prof_options.solvers = solvers;
prof_options.outdir = outdir;
prof_options.stamp = stamp;
prof_options.feature_and_time = feature_and_time;

prec = (1:12);
nprec = length(prec);
tau = 10.^(-prec);
prof_output = cell(1, 2*nprec);
format long
for iprec = 1 : 2*nprec
    % If `natural_stop` is true, the number of function evaluations is the amount used by the solver
    % when it stops naturally.
    prof_options.natural_stop = (iprec > nprec);
    real_iprec = mod(iprec-1, nprec) + 1;  % The real index of the precision.
    prof_options.tau = tau(real_iprec);
    prof_output{iprec} = perfprof(frec, fmin, prof_options);
    %dataprof(frec, fmin, pdim, prof_options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appearance of the plots.
%fontsize = 12;
linewidth = 1;
bleu=[0.16, 0.4470, 0.7410];
vert=[0, 0.6, 0];
colors  = {bleu, 'k', 'b', 'r', vert, bleu, 'k', 'b', 'r', vert};
%lines   = {'-', '-.', '--', ':', '-', '-.', '--', ':', '-', '-.'};
lines   = {'-', '-', '-', '-', '-', '-', '-', '-', '-', '-'};
hfig = figure("visible", false, 'DefaultAxesPosition', [0, 0, 1, 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iprec = 1 : 2*nprec
    subplot(2, nprec, iprec);
    if iprec > nprec
        lw = linewidth;
    else
        lw = linewidth * 3;
    end
    for is = 1 : ns
       plot(prof_output{iprec}.profile{is}(1,:), prof_output{iprec}.profile{is}(2,:), ...
           lines{is}, 'Color', colors{is},  'Linewidth', lw);
       hold on;
    end
    xlabel(sprintf('%d', iprec));
    axis([0 prof_output{iprec}.cut_ratio 0 1]);
    grid on;
    %pbaspect([1 1 1]);
end
for iprec = 1 : 2*nprec
    ha = get(gcf,'children');
    real_iprec = mod(iprec-1, nprec) + 1;  % The real index of the precision.
    %if iprec > nprec
    %    set(ha(iprec),'position', [0.01+(nprec-real_iprec)/(nprec), 0.53, 0.9/(nprec), 0.4]);
    %else
    %    set(ha(iprec),'position', [0.01+(nprec-real_iprec)/(nprec), 0.1, 0.9/(nprec), 0.4]);
    %end
    if iprec > nprec
        set(ha(iprec),'position', [0.01+(nprec-real_iprec)/(nprec), 0.53, 0.8/(nprec), 0.35]);
    else
        set(ha(iprec),'position', [0.01+(nprec-real_iprec)/(nprec), 0.1, 0.8/(nprec), 0.35]);
    end
end

% The following appears only in the last subplot, but it is sufficient for our use.
ylabel(strrep(test_feature, '_', '\\_'));
legend(solvers,'Location', 'southeast','Orientation','vertical');

% Save the figure as eps.
figname = strcat(stamp, '.', 'summary', '.', feature_and_time);
epsname = fullfile(outdir, strcat(figname,'.eps'));
set(gcf,'position',[0, 0, 4600, 920]);
saveas(hfig, epsname, 'epsc2');

% Try converting the eps to pdf.
try
    system(['epstopdf ',epsname]);
    summary_file = fullfile(outdir, strcat(figname,'.pdf'));
catch
    summary_file = epsname;
end
fprintf('\nSummary for problem type %s with test feature %s:\n\n%s\n\n', options.type, test_feature, summary_file);


% For convenience, save a copy of `problems.txt` and the figures in data_dir. They will be
% replaced in next test with the same `solvers` and `dimrange`.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the old files.
delete(fullfile(data_dir, strcat(stamp, '.*.problems.txt')));
delete(fullfile(data_dir, strcat(stamp, '.perf_*.pdf')));
delete(fullfile(data_dir, strcat(stamp, '.perf_*.eps')));
delete(fullfile(data_dir, strcat(stamp, '.summary.*.eps')));
delete(fullfile(data_dir, strcat(stamp, '.summary.*.pdf')));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
copyfile(fprob, data_dir);
epsfiles = dir(fullfile(outdir, '*.eps'));
for k = 1 : length(epsfiles)
    epsname = epsfiles(k).name;
    pdfname = strrep(epsname, '.eps', '.pdf');
    if exist(fullfile(outdir, pdfname), 'file')
        copyfile(fullfile(outdir, pdfname), data_dir);
    else
        copyfile(fullfile(outdir, epsname), data_dir);
    end
end
