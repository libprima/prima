function info = getMexLibgcc()
%GETMEXLIBGCC Resolve libgcc_s used by the current configuration of Fortran MEX and its GCC version.
% Output struct fields:
%   libgccPath: full path to the resolved libgcc_s.so.1; empty if not detectable
%   libgccStrings: strings embedded in libgcc_s.so.1 (a string separated by line breaks); empty if not detectable
%   latestGccVersion: latest GCC version string embedded in libgcc_s.so.1 (e.g., '14.0.0'); empty if not detectable
%
% Requirements: Linux, ldd, strings, grep, tail. No root needed.

    % Get name of the current function
    callstack = dbstack;
    funName = callstack(1).name;

    % Preconditions
    assert(isunix && ~ismac, sprintf('%s: This function targets Linux.', funName));
    assert(isCommandAvailable('ldd'), sprintf('%s: ldd command not found.', funName));
    assert(isCommandAvailable('strings'), sprintf('%s: strings command not found.', funName));
    assert(isCommandAvailable('grep'), sprintf('%s: grep command not found.', funName));
    assert(isCommandAvailable('tail'), sprintf('%s: tail command not found.', funName));

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

    libgccStrings = '';
    latestGccVersion = '';
    if ~isempty(libgccPath)
        % Extract strings from libgcc
        [~, libgccStrings] = system(sprintf('strings "%s"', libgccPath));

        % Find the latest GCC version mentioned in the strings, supposing that the latest is mentioned in the last line
        [~, latestGccString] = system(sprintf('strings "%s" | grep -E "^GCC_[0-9]+" | tail -n1', libgccPath));
        latestGccString = strtrim(latestGccString);
        % Try to parse a version like 14.0.0 from the latestGccString
        tok = regexp(latestGccString, '([0-9]+(\.[0-9]+)?+(\.[0-9]+)?)', 'tokens', 'once');
        if ~isempty(tok)
            latestGccVersion = tok{1};
        end
    end

    info = struct( ...
        'libgccPath', libgccPath, ...
        'libgccStrings', libgccStrings, ...
        'latestGccVersion', latestGccVersion);
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
