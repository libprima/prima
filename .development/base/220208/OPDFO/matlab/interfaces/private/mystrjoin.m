function str = mystrjoin(C, delimiter) 
% MYSTRJOIN replaces strjoin, which was introduced only after R2013a.
%
%   ***********************************************************************
%   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk) 
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%
%   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: private (not supposed to be called by users)
%
% Remarks: None
%
% TODO: None
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mystrjoin starts
callstack = dbstack; 
funname = callstack(1).name; % Name of the current function

if nargin < 1 
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: at least 1 input.', funname);
elseif nargin == 1
    delimiter = ' ';
end

if ~isa(C, 'cell') || ~isvector(C) || any(~cellfun(@ischarstr, C))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), ...
    '%s: UNEXPECTED ERROR: input 1 should be a cell array of character vectors or strings.', funname);
end
    
len = length(C);
if (len == 0)
    str = '';
elseif (len == 1)
    str = C{1};
else
    CD = cell(2*length(C)-1, 1);
    CD(2*(1:len)-1) = C;
    CD(2*(1:len-1)) = {delimiter};
    str = sprintf('%s', CD{:});
end

% mystrjoin ends
return

function ics = ischarstr(x) 
ics = isa(x,'char') || isa(x, 'string') || isempty(x);
return
