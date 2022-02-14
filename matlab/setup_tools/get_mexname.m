function mexname = get_mexname(solver, precision, debug_flag, variant)
%GET_MEXNAME returns the name of the mexified `solver` according to `precision`, `debug_flag`, and
% `variant`.

switch solver
case {'gethuge'}
    mexname = [solver, '_', precision(1)];
otherwise
    mexname = [solver, '_', precision(1), dbgstr(debug_flag), variant(1)];
end

return


function s = dbgstr(debug_flag)
%DBGSTR returns a string according to `debug_flag`: 'g' for true and 'O' for false.
if (debug_flag)
    s = 'g';
else
    s = 'O';
end
return
