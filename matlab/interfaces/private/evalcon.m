function [cineq, ceq, succ] = evalcon(invoker, nonlcon, x)
%EVALCON evaluates a constraint function `[cineq, ceq] = nonlcon(x)`.
% In particular, it uses a 'moderated extreme barrier' to cope with 'hidden constraints' or
% evaluation failures.

succ = true;

try
    [cineq, ceq] = nonlcon(x);
catch
    succ = false;
    cineq = NaN;
    ceq = NaN;
end

if ~(isempty(cineq) || isnumeric(cineq))
    succ = false;
    cineq = NaN;
end

if ~(isempty(ceq) || isnumeric(ceq))
    succ = false;
    ceq = NaN;
end

% Use a 'moderated extreme barrier' to cope with 'hidden constraints' or evaluation failures.
hugecon = gethuge('con');

if any(isnan(cineq) | ~isreal(cineq) | cineq > hugecon)
    wid = sprintf('%s:ConstraintFailure', invoker);
    xstr = sprintf('%f  ', x);
    if any(~isreal(cineq))
        cstr = sprintf('%f%+fi  ', [real(cineq(:)), imag(cineq(:))].');
    else
        cstr = sprintf('%f  ', cineq);
    end
    wmsg = sprintf('%s: Constraint function returns cineq =\n%s\nAny value that is not real or above hugecon = %1.2e is replaced by hugecon.\nThe value of x is:\n%s\n', invoker, cstr, hugecon, xstr);
    warning(wid, '%s', wmsg);
    %warnings = [warnings, wmsg];  % We do not record this warning in the output.

    % Apply the moderated extreme barrier:
    cineq(~isreal(cineq) | cineq~= cineq | cineq > hugecon) = hugecon;
end

if any(isnan(ceq) | ~isreal(ceq) | abs(ceq) > hugecon)
    wid = sprintf('%s:ConstraintFailure', invoker);
    xstr = sprintf('%f  ', x);
    if any(~isreal(ceq))
        cstr = sprintf('%f%+fi  ', [real(ceq(:)), imag(ceq(:))].');
    else
        cstr = sprintf('%f  ', ceq);
    end
    wmsg = sprintf('%s: Constraint function returns ceq =\n%s\nAny value that is not real or with an absolute value above hugecon = %1.2e is replaced by hugecon.\nThe value of x is:\n%s\n', invoker, cstr, hugecon, xstr);
    warning(wid, '%s', wmsg);
    %warnings = [warnings, wmsg];  % We do not record this warning in the output.

    % Apply the moderated extreme barrier:
    ceq(~isreal(ceq) | ceq ~= ceq | ceq > hugecon) = hugecon;
    ceq(ceq < -hugecon) = -hugecon;
end

% This part is NOT an extreme barrier. We replace extremely negative values of
% cineq (which leads to no constraint violation) by -hugecon. Otherwise,
% NaN or Inf may occur in the interpolation models.
cineq(cineq < -hugecon) = -hugecon;

cineq = double(real(cineq(:))); % Some functions like 'asin' can return complex values even when it is not intended
ceq = double(real(ceq(:)));
return
