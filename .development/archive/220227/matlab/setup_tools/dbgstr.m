function s = dbgstr(debug_flag)
%DBGSTR returns a string according to `debug_flag`: 'g' for true (debug) and 'O' for false (not debug).
% Any strings that can distinguish true and false will work. We choose 'g' and 'O' because MEX uses
% `-g` as the flag for debugging and `-O` for no debugging but optimizing.

if (debug_flag)
    s = 'g';
else
    s = 'O';
end

return
