function T = rmempty(S) % Remove empty fields in a structure
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if isempty(S)
    S = struct(); % Here we do not distinguish empty objects. It is fine in this package, but may not be in others
end
if ~isa(S, 'struct')
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: input should be a structure.', funname);
end
fn = fieldnames(S);
empty_index = cellfun(@(f) isempty(S.(f)), fn);
T = rmfield(S, fn(empty_index));
return
