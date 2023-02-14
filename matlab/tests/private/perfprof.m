function output = perfprof(frec, fmin, options)
%This function plots the performance profiles of solvers.
% frec: trajectory of function values; frec(ip, is, ir, k) is the function value of the ip-th problem
% obtained by the is-th solver at the ir-th random run at the k-th iteration.
% fmin: the minimal function values; fmin(ip) is the minimal function value of the ip-th problem.
% tau: the tolerance of convergence.
% solvers: the list of solvers.

% Parameters.
penalty = 100;
cut = 1.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appearance of the plots.
fontsize = 12;
linewidth = 1;
% Colors.
bleu=[0.16, 0.4470, 0.7410];
vert=[0, 0.6, 0];
colors  = {bleu, 'k', 'b', 'r', vert, bleu, 'k', 'b', 'r', vert};
lines   = {'-', '-.', '--', ':', '-', '-.', '--', ':', '-', '-.'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[np, ns, nr, maxfun] = size(frec);
M = maxfun;

% T(ip, is, ir) is the number of function evaluations that the is-th solver needs to solve the ip-th
% problem (up to tolerance tau) at the ir-th random run.
T = NaN(np, ns, nr);

f0 = -Inf(np, nr);
for ip = 1:np
    for ir = 1:nr
        f0(ip,ir) = frec(ip, 1, ir, 1);
    end
end

tau = options.tau;
for ip = 1:np
    for is = 1:ns
        for ir = 1:nr
            fthreshold = tau*f0(ip,ir) + (1-tau)*fmin(ip);
            % We need to ensure fthreshold >= fmin(ip), which may not be true due to rounding
            % errors when fmin(ip)=f0(ip,ir).
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fthreshold = max(fthreshold, fmin(ip));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (min(frec(ip, is, ir, 1:M)) <= fthreshold)
                % Do not change the "if .. else ..." order, as frec(ip, is, ir, 1:M) may be all NaNs.
                T(ip, is, ir) = find(frec(ip, is, ir, 1:M) <= fthreshold, 1, 'first');
            else
                T(ip, is, ir) = NaN;
            end
        end
    end
end

pp = cell(ns, nr);
cut_ratio = 0;
for ir = 1 : nr
    Tmin = min(T(:, :, ir), [], 2);

    % Calculate r. r(ip, is, ir) is the number of function evaluations that the is-th solvers needs
    % to solve the ip-th problem during the ir-th run, normalized by the minimal number of function
    % evaluations needed by all the solvers for this problem during the same run. We replace r by
    % log2(r) for better scaling. Why not using semilogx instead? Because semilogx only supports log10!
    r = zeros(np, ns);
    for ip = 1 : np
      r(ip, :) = T(ip, :, ir)/Tmin(ip);
    end
    r = log2(r);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    max_ratio = max(1.0D-1, max(max(r)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    penalty_ratio = penalty*max_ratio;
    cut_ratio = max(cut_ratio, cut*max_ratio);
    r(isnan(r)) = penalty_ratio;
    r = sort(r);

    for is = 1:ns
        [xx, yy] = stairs(r(:, is), (1:np)/np);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The following ensures the correctness of the profiles in extreme
        % cases like one of the solvers always performs the best, or one of
        % the solvers cannot solve any problem.
        %xx = [0; xx(1); xx; penalty_ratio];
        %y = [0; 0; yy; yy(end)];
        xx = [0; xx(1); xx];
        yy = [0; 0; yy];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xx = xx + 10*max(xx)*eps*(0: length(xx)-1)';
        pp{is, ir} = [xx'; yy'];
    end
end


x = cell(1, ns);
for is = 1 : ns
    x{is} = [];
    for ir = 1 : nr
        x{is} = [x{is}, pp{is, ir}(1, :)];
    end
    x{is} = sort(x{is});
end

y = cell(1, ns);
perf_prof = cell(1, ns);
for is = 1 : ns
    r = NaN(nr, length(x{is}));
    for ir = 1 : nr
        for ix = 1 : length(x{is})
            r(ir, ix) = pp{is, ir}(2, find(pp{is, ir}(1, :) <= x{is}(ix), 1, 'last'));
        end
    end
    y{is} = mean(r, 1);
    perf_prof{is} = [x{is}; y{is}];
end

clf;
hfig=figure("visible", false);  % Plot the figure without displaying it.
for is = 1:ns
    plot(perf_prof{is}(1, :), perf_prof{is}(2, :), lines{is}, 'Color', colors{is},  'Linewidth', linewidth);
    hold on;
end

%output.success_rate = success_rate;
output.profile = perf_prof;
output.cut_ratio = cut_ratio;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([0 cut_ratio 0 1]);
yticks(0 : 0.1 : 1);
yticklabels({'0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Legends and title should be added.
%title(sprintf('Performance Profile with tolerance $\\tau=10^{%d}$', int32(log10(tau))), 'interpreter', 'latex');
solvers = options.solvers;
for is = 1:ns
    solvers{is} = regexprep(solvers{is}, '_4test', '');
    solvers{is} = regexprep(solvers{is}, '_classical$', ' (classical)');
    solvers{is} = regexprep(solvers{is}, '_single$', ' (single)');
    solvers{is} = regexprep(solvers{is}, '_quadruple$', ' (quadruple)');
    %solvers{is} = regexprep(solvers{is}, 'newuoa', 'NEWUOA');
end
if (ns > 3)
    legend(solvers,'Location', 'southeast','Orientation','vertical');
else
    %legend(solvers,'Location', 'northoutside','Orientation','horizontal');
    legend(solvers,'Location', 'southeast','Orientation','vertical');
end

xlabel('$\log_2(\alpha), \quad \alpha = \mathrm{NF}/\mathrm{NF}_{\min}$', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('$\pi_s(\alpha)$', 'fontsize', fontsize, 'interpreter', 'latex');
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .0175 0])
set(gca,'FontSize',fontsize);

% Save the figure as eps.
figname = strcat(options.stamp, '.', 'perf_', int2str(int32(-log10(tau))), '.', options.feature_and_time);
epsname = fullfile(options.outdir, strcat(figname,'.eps'));
saveas(hfig, epsname, 'epsc2');

% Try converting the eps to pdf.
try
    system(['epstopdf ',epsname]);
catch
    % Do nothing in case of failure.
end
