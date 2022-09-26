function [f, succ] = evalobj(invoker, fun, x, hugefun)
%EVALOBJ evaluates an objective function `f = fun(x)`.
% In particular, it uses a 'moderated extreme barrier' to cope with 'hidden constraints' or
% evaluation failures.

succ = true;

try
    f = fun(x);
catch exception
    succ = false;
    f = NaN;
    wid = sprintf('%s:ObjectiveFailure', invoker);
    wmsg = sprintf('%s: Objective function fails with the following error:\n  %s: %s\n  Error occurred in %s, line %d', ...
        invoker, exception.identifier, exception.message, exception.stack(1).file, exception.stack(1).line);
    warning(wid, '%s', wmsg);
end

if numel(f) ~= 1
    % succ = false;  % Not needed if we raise an error.
    % Public/normal error
    error(sprintf('%s:ObjectiveNotScalar', invoker), '%s: objective function should return a scalar value.', invoker);
end

if ~isnumeric(f)
    succ = false;
    f = NaN;
end

% Use a 'moderated extreme barrier' to cope with 'hidden constraints' or evaluation failures.
if isnan(f) || ~isreal(f) || f > hugefun
    wid = sprintf('%s:ObjectiveAbnormalReturn', invoker);
    xstr = sprintf('%g    ', x);
    if ~isreal(f)
        fstr = sprintf('%g%+gi', real(f), imag(f));
    else
        fstr = sprintf('%g', f);
    end
    wmsg = sprintf('%s: Objective function returns %s, which is replaced by hugefun = %g.\nThe value of x is:\n%s\n', invoker, fstr, hugefun, xstr);
    warning(wid, '%s', wmsg);
    %warnings = [warnings, wmsg];  % We do not record this warning in the output.

    % Apply the moderated extreme barrier:
    f = hugefun;
end

f = double(real(f)); % Some functions like 'asin' can return complex values even when it is not intended
return
