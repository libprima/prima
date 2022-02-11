function [f, succ] = evalobj(invoker, fun, x)
%EVALOBJ evaluates an objective function `f = fun(x)`.
% In particular, it uses a 'moderated extreme barrier' to cope with 'hidden constraints' or
% evaluation failures.

succ = true;

try
    f = fun(x);
catch
    succ = false;
    f = NaN;
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
hugefun = gethuge('fun');
if isnan(f) || ~isreal(f) || f > hugefun
    wid = sprintf('%s:ObjectiveFailure', invoker);
    xstr = sprintf('%f  ', x);
    if ~isreal(f)
        fstr = sprintf('%f%+fi', real(f), imag(f));
    else
        fstr = sprintf('%f', f);
    end
    wmsg = sprintf('%s: Objective function returns %s, which is replaced by hugefun = %1.2e.\nThe value of x is:\n%s\n', invoker, fstr, hugefun, xstr);
    warning(wid, '%s', wmsg);
    %warnings = [warnings, wmsg];  % We do not record this warning in the output.

    % Apply the moderated extreme barrier:
    f = hugefun;
end

f = double(real(f)); % Some functions like 'asin' can return complex values even when it is not intended
return
