function hugenum = gethuge(data_type, precision)

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

if nargin < 2
    precision = 'double';
end

if ~ischarstr(precision) || ~ismember(precision, all_precisions())
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid precision received', funname);
end

gethuge_prec = str2func(get_mexname('gethuge', precision));

hugenum = gethuge_prec(data_type);

return
