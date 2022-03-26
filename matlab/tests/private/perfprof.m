function output = perfprof(frec, fmin, options)
%This function plots the performance profiles of solvers.
% frec: trajectory of function values; frec(ip, is, ir, k) is the function value of the ip-th problem
% obtained by the is-th solver at the ir-th random run at the k-th iteration.
% fmin: the minimal function values; fmin(ip) is the minimal function value of the ip-th problem.
% tau: the tolerance of convergence.
% solvers: the list of solvers.

% Parameters.
delsame = 0;
penalty = 2;
cut = 1.05;
fontsize = 12;
linewidth = 1;
% Colors.
bleu=[0.16, 0.4470, 0.7410];
vert=[0, 0.6, 0];
colors  = {bleu, 'k', 'b', 'r', vert, bleu, 'k', 'b', 'r', vert};
lines   = {'-', '-.', '--', ':', '-', '-.', '--', ':', '-', '-.'};

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

% T(ip, is) is the average number of function evaluations that the is-solver needs to solve the
% ip-th problem (up to the tolerance tau), average taken across the random runs.
T = mean(T, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(delsame==1)
warning('Deleting all the problems for which all the solvers perform the same.');
% Delete the problems for which all the solvers performs the same.
    mask = true(np,1);
    for ip = 1:np
        if (sum(isnan(T(ip,:))) == ns || (sum(isnan(T(ip,:))) == 0 && max(T(ip,:)) <= min(T(ip,:))+1))
            mask(ip) = false;
        end
    end
    T = T(mask,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[np, ns] = size(T);
if(np == 0)  % Prevent T from being empty.
    T=ones(1,ns);
    np =1;
end
Tmin = min(T, [], 2);

% Calculate r. r(ip, is) is the number of function evaluations that the is-th solvers needs to solve
% the ip-th problem, normalized by the minimal number of function evaluations needed by all the
% solvers for this problem. We replace r by log2(r) for better scaling. Why not using semilogx
% instead? Because semilogx only supports log10!
r = zeros(np, ns);
for ip = 1: np
  r(ip, :) = T(ip, :)/Tmin(ip);
end
r = log2(r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_ratio = max(1.0D-1,max(max(r)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
penalty_ratio = penalty*max_ratio;
cut_ratio = cut*max_ratio;
r(isnan(r)) = penalty_ratio;
r = sort(r);

success_rate = cell(1, ns);
perf_prof = cell(1, ns);

clf;
hfig=figure("visible", false);  % Plot the figure without displaying it.
for is = 1:ns
    [x,y] = stairs(r(:,is), (1:np)/np);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following ensures the correctness of the profiles in extreme
    % cases like one the solvers alway performs the best, or one of
    % the solvers cannot solve any problem.
    x = [0; x(1); x; penalty_ratio];
    y = [0; 0; y; y(end)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    perf_prof{is} = [x'; y'];
    plot(x, y, lines{is}, 'Color', colors{is},  'Linewidth', linewidth);
    hold on;

    % Evaluate the success rate.
    if ~any(x > 0 & x < penalty_ratio)
        success_rate{is} = ones(1,3);
        continue;
    end
    xx = x(x>0 & x < penalty_ratio);
    yy = y(x>0 & x < penalty_ratio);
    xx = [0; xx; xx(end)*cut];
    yy = [yy(1); yy; yy(end)];
    success_rate{is} = NaN(1,3);
    success_rate{is}(1) = yy(1);  % The success rate corresponding to NF/NFMIN = 1 (i.e., NF=NFMIN).
    success_rate{is}(3) = yy(end);  % The success rate with "infinite budget".


    % The following lines calculates the average success rate.
    k = floor(length(xx)/2);  % In theory, length(xx) is even
    % The following line calculates the simple average
    %success_rate(is, 2) = sum(yy(2*(1:k)-1).*(xx(2*(1:k))-xx(2*(1:k)-1)))/xx(end);
    % The following line calculates the weighted average with the weight = exp(-NF/NFMIN), putting
    % more weight when NF/NFMIN is small.
    success_rate{is}(2) = sum(yy(2*(1:k)-1).*(xx(2*(1:k))-xx(2*(1:k)-1)).*exp(-xx(2*(1:k)-1)))/sum((xx(2*(1:k))-xx(2*(1:k)-1)).*exp(-xx(2*(1:k)-1)));
end

output.success_rate = success_rate;
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
figname = strcat(options.stamp, '.', 'perf_', int2str(int32(-log10(tau))), '.', options.time);
epsname = fullfile(options.outdir, strcat(figname,'.eps'));
saveas(hfig, epsname, 'epsc2');

% Try converting the eps to pdf.
try
    system(['epstopdf ',epsname]);
catch
    % Do nothing in case of failure.
end
