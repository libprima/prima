function info = getMexLibgcc()
%GETMEXLIBGCC Resolve libgcc_s used by the current configuration of Fortran MEX and its GCC version.
% Output struct fields:
%   libgccPath      - full path to the resolved libgcc_s.so.1; empty if not detectable
%   gccVersion      - GCC version string (e.g., '14.2.0'); empty if not detectable
%
% Requirements: Linux, ldd, strings, grep, head. No root needed.

    % Get name of the current function
    callstack = dbstack;
    funName = callstack(1).name;

    % Preconditions
    assert(isunix && ~ismac, sprintf('%s: This function targets Linux.', funName));
    assert(isCommandAvailable('ldd'), sprintf('%s: ldd command not found.', funName));
    assert(isCommandAvailable('strings'), sprintf('%s: strings command not found.', funName));
    assert(isCommandAvailable('grep'), sprintf('%s: grep command not found.', funName));
    assert(isCommandAvailable('head'), sprintf('%s: head command not found.', funName));

    % exampleFile is an example provided by MATLAB for trying MEX.
    % NOTE: MATLAB MAY CHANGE THE LOCATION OF THIS FILE IN THE FUTURE.
    exampleFile = fullfile(matlabroot, 'extern', 'examples', 'refbook', 'timestwo.F');
    [~, mexName] = fileparts(exampleFile);
    if ~exist(exampleFile, 'file')
        error('%s: The official MEX example not found: %s', funName, exampleFile);
    end

    % Build MEX into a temp directory
    % NOTE: We do not try setting up MEX here, because we want to probe the current MEX setup.
    outDir = tempname();
    mkdir(outDir);
    c = onCleanup(@() safeCleanup(outDir));  % Cleanup temp dir on function exit
    try
        clear('timestwo');
        evalc('mex(''-outdir'', outDir, ''-output'', mexName, exampleFile)');
    catch ME
        error('%s: Failed to build MEX from %s.\nThe error message is\n%s\nMake sure that MEX is properly set up.', funName, exampleFile, ME.message);
    end

    % Find produced MEX
    mexFile = fullfile(outDir, [mexName, '.', mexext()]);
    if ~exist(mexFile, 'file')
        error('%s: Could not find built MEX: %s', funName, mexFile);
    end

    % Resolve which libgcc_s this MEX will use
    [st, lddOut] = system(sprintf('ldd "%s"', mexFile));
    if st ~= 0
        error('%s, ldd failed: %s', funName, lddOut);
    end
    libgccPath = parseLibPath(lddOut, 'libgcc_s.so.1');

    % Extract GCC version string
    gccVersion = '';
    if ~isempty(libgccPath)
        % We suppose that the system embeds a GCC build marker the strings tool can find.
        [~, marker] = system(sprintf('strings "%s" | grep -Ei "^GCC:|GCC \\(" | head -n1', libgccPath));
        marker = strtrim(marker);
        % Try to parse a version like 14.2.0 from the marker
        % Examples:
        %   GCC: (Ubuntu 14.2.0-3ubuntu1) 14.2.0
        %   GCC: (GNU) 13.2.1 20240109
        %   GCC: (Red Hat 13.2.1-7) 13.2.1 20240316
        tok = regexp(marker, '([0-9]+\.[0-9]+(\.[0-9]+)?)', 'tokens', 'once');
        if ~isempty(tok)
            gccVersion = tok{1};
        end
    end

    info = struct( ...
        'libgccPath', libgccPath, ...
        'gccVersion', gccVersion);
end


function p = parseLibPath(lddOut, soname)
    p = '';
    lines = regexp(lddOut, '\r?\n', 'split');
    for i = 1:numel(lines)
        L = strtrim(lines{i});
        if startsWith(L, soname)
            % Format: libgcc_s.so.1 => /path/to/libgcc_s.so.1 (0x....)
            tok = regexp(L, '=>\s+(\S+)\s+\(', 'tokens', 'once');
            if ~isempty(tok)
                p = tok{1};
            end
            break;
        end
    end
end

function safeCleanup(p)
    if ~exist(p, 'dir')
        return;
    end
    try
        rmdir(p, 's');
    catch
    end
end

function isAvailable = isCommandAvailable(cmd)
    [status, location] = system(['command -v ', cmd]);
    isAvailable = (status == 0 && ~isempty(strtrim(location)));
end
