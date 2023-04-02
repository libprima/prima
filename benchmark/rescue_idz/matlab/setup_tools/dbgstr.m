function s = dbgstr(debug_flag)
%DBGSTR returns a string according to `debug_flag`: 'g' for true (debug) and 'O' for false (not debug).
% This function is used at two places:
% - compile.m, where `dbgstr` is used to generate the MEX options (MEX uses `-g` as the flag for
%   debugging and `-O` for no debugging but optimizing);
% - get_mexname.m, where `dbgstr` is used to generate names of the MEX functions.

if (debug_flag)
    s = 'g';
else
    s = 'O';
end

return
