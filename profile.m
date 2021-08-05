function T = profile(frec, fmin, tau, n, testfeature)
%
% This version is not inteded to be released. It is only for test. 
%
% All rights reserved. 
%
% ZHANG Zaikun, 08/08/2016
% Department of Applied Mathematics, The Hong Kong Polytechnic University

delsame = 0;
penalty = 2;
cut = 1.1;
fontsize = 12;
linewidth = 1;

bleu=[0.16, 0.4470, 0.7410];
vert=[0, 0.6, 0];
colors  = {bleu, 'k', 'b', 'r', vert, bleu, 'k', 'b', 'r', vert};
lines   = {'-', '-.', '--', ':', '-', '-.', '--', ':', '-', '-.'};

[np, ns, nr, maxfun] = size(frec);
M = maxfun;

T = NaN(np, ns, nr);

f0 = -Inf(np, nr);
for ip = 1:np
    for ir = 1:nr 
        f0(ip,ir) = frec(ip, 1, ir, 1);
    end
end

for ip = 1:np
    for is = 1:ns
        for ir = 1:nr
            fthreshold = tau*f0(ip,ir) + (1-tau)*fmin(ip);
            % We need to ensure fthreshold >= fmin(ip), which may not be true 
            % due to rounding errors when fmin(ip)=f0(ip,ir). 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fthreshold = max(fthreshold, fmin(ip)); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (min(frec(ip, is, ir, 1:M)) <= fthreshold) % Do not change the "if .. else ..." order, because frec(ip, is, ir, 1:M) may be a vector of NaNs.
                T(ip, is, ir) = find(frec(ip, is, ir, 1:M) <= fthreshold, 1, 'first');
            else
                T(ip, is, ir) = NaN;
            end
        end
    end
end

T = mean(T, 3);

sol = textread('solvers', '%s');
for is = 1:ns
    sol{is} = regexprep(sol{is}, '_4test', '');
    %sol{is} = regexprep(sol{is}, 'newuoa', 'NEWUOA');
end

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

clf;
hfig=figure(1);
for is = 1:ns
    [xs,ys] = stairs(r(:,is), (1:np)/np);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following ensures the correctness of the profiles in extreme
    % cases like one the solvers always performs the best, or one of
    % the solvers cannot solve any problem.
    xs = [0; xs(1); xs; penalty_ratio];
    ys = [0; 0; ys; ys(end)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(xs, ys, lines{is}, 'Color', colors{is},  'Linewidth',linewidth);
    hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([0 cut_ratio 0 1]);
yticks(0 : 0.1 : 1);
yticklabels({'0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Legends and title should be added.
%title(sprintf('Performance Profile with tolerance $\\tau=10^{%d}$', int32(log10(tau))), 'interpreter', 'latex');
if (ns > 3)
    legend(sol,'Location', 'southeast','Orientation','vertical');
else
    %legend(sol,'Location', 'northoutside','Orientation','horizontal');
    legend(sol,'Location', 'southeast','Orientation','vertical');
end

xlabel('$\log_2(\alpha), \quad \alpha = \mathrm{NF}/\mathrm{NF}_{\min}$', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('$\pi_s(\alpha)$', 'fontsize', fontsize, 'interpreter', 'latex');
xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .0175 0])
set(gca,'FontSize',fontsize);

set(hfig,'visible','off');
figname = strcat(testfeature, '_', int2str(n), '_', int2str(int32(-log10(tau))),'_perf');
epsname = strcat(figname,'.eps');
saveas(hfig, epsname, 'epsc2');
system('mkdir -p results');
system(['epstopdf ',epsname]);
system(['mv ', epsname, ' ./results']);
system(['mv ', strcat(figname, '.pdf'), ' ./results']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = NaN(ns, M);
for is = 1:ns
    for i = 1:M
        D(is, i) = length(find(T(:, is) <= i))/np;
    end
    for i = M+1:1.2*M
        D(is, i) = D(is, M);
    end
end

clf;
hfig=figure(2);
for is = 1:ns
    [xs, ys] = stairs([1:1.2*M]/(n+1), D(is, :));
    plot(xs, ys, lines{is}, 'Color', colors{is},  'Linewidth',linewidth);
    hold on;
end
axis([ 0 1.1*M/(n+1) 0 1 ]);
yticks(0 : 0.1 : 1);
yticklabels({'0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'});

% Legends and title should be added.
if (ns > 3)
    legend(sol,'Location', 'eastoutside','Orientation','vertical');
else
    %legend(sol,'Location', 'northoutside','Orientation','horizontal');
    legend(sol,'Location', 'southeast','Orientation','vertical');
end

xlabel('$\beta = \mathrm{NF}/(n+1)$', 'fontsize',fontsize, 'interpreter', 'latex');
ylabel('$\delta_s(\beta)$', 'fontsize', fontsize, 'interpreter', 'latex');
xlabh = get(gca, 'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .0175 0])
set(gca, 'FontSize', fontsize);

set(hfig,'visible','off');
figname = strcat(testfeature, '_', int2str(n), '_', int2str(int32(-log10(tau))),'_data');
epsname = strcat(figname,'.eps');
saveas(hfig, epsname, 'epsc2');
system('mkdir -p results');
system(['epstopdf ',epsname]);
system(['mv ', epsname, ' ./results']);
system(['mv ', strcat(figname, '.pdf'), ' ./results']);
