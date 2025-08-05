function information = package_info(request)
%PACKAGE_INFO returns information about the package.
%
%   ***********************************************************************
%   Author:     Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%
%   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: private (not supposed to be called by users)
%
% TODO: None
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% package_info starts

% Who is calling this function? Is it a correct invoker?
invoker_list = {'prima_norma', 'preprima_norma'};
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if (length(callstack) == 1 || ~ismember(callstack(2).name, invoker_list))
    % Private/unexpected error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s should only be called by %s.', funname, funname, strjoin(invoker_list, ', '));
else
    invoker = callstack(2).name; % Name of the function who calls this function
end

name = 'PRIMA';

about = 'PRIMA: Reference Implementation for Powell''s methods with Modernization and Amelioration. PRIMA provides the reference implementation of late Professor M. J. D. Powell''s derivative-free optimization methods, namely COBYLA, UOBYQA, NEWUOA, BOBYQA, and LINCOA.';

author = 'Zaikun Zhang';

email = 'zaikun.zhang@polyu.edu.hk';

url = 'www.libprima_norma.net';

maintainer = 'Zaikun Zhang';

credits = {'Tom M. Ragonneau', 'Zaikun Zhang'};

license = 'LGPLv3+';

version = '0.9';

date = 'January, 2023';

status = 'Development';

message = 'Dedicated to late Professor M. J. D. Powell FRS (29 July 1936--19 April 2015).';

copyright = sprintf('Copyright 2020--%d, Zaikun Zhang', year(datetime()));

switch lower(request)
case 'name'
    information = name;
case 'about'
    information = about;
case 'author'
    information = author;
case 'email'
    information = email;
case 'url'
    information = url;
case 'maintainer'
    information = maintainer;
case 'credits'
    information = credits;
case 'copyright'
    information = copyright;
case 'license'
    information = license;
case 'version'
    information = version;
case 'date'
    information = date;
case 'status'
    information = status;
case 'message'
    information = message;
case {'info', 'information'}
    information = struct('name', name, 'about', about, 'author', author, 'email', email, 'url', url, 'maintainer', maintainer, 'credits', [], 'copyright', copyright, 'license', license, 'version', version, 'date', date, 'status', status, 'message', message);
    % information = struct(..., 'credits', credits, ...) will produce
    % a cell array of size 1x2, which is not desired. Thus we fist
    % define information with information.credits = [], and then assign
    % the following value:
    information.credits = credits;
otherwise % Public/expected error
    error(sprintf('%s:UnrecognizedString', invoker), '%s: unrecognized string received.', invoker);
end

% package_info ends
return
