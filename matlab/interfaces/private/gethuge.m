function hugenum = gethuge(data_type, precision)
%GETHUGE calls returns a huge number according to `data_type` and `precision`.

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

if nargin < 1 || nargin > 2
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid number of inputs.', funname);
else
    if ~ischarstr(data_type)
        error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid data_type received.', funname);
    end
    if nargin == 2 && ~(ischarstr(precision) && ismember(precision, all_precisions()))
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid precision received.', funname);
    end
end

if nargin < 2
    precision = 'double';
end

gethuge_prec = str2func(get_mexname('gethuge', precision));

hugenum = gethuge_prec(data_type);

return
