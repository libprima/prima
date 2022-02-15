function freeform(files)
%FREEFORM refactors the Fortran code in `files` from the fixed form to the free form.
%
% Coded by Zaikun ZHANG (www.zhangzk.net) in August, 2020.

if nargin < 1 || isempty(files) || strcmpi(files, 'ALL')
    listing = dir();
    files = {listing.name};
else
    files = {files};
end

for ifile = 1 : length(files)
    filename = files{ifile};
    if endsWith(filename, '.f') || endsWith(filename, '.f90')
        fid = fopen(filename, 'r');  % Open file for reading.
        if fid == -1
            error('Cannot open file %s', filename);
        end
        data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
        fclose(fid);
        cstr = data{1};
        cstr(count(cstr, 'C') == strlength(cstr)) = replace(cstr(count(cstr, 'C') == strlength(cstr)), 'C', '!');
        for jc = 1 : length(cstr)
            if startsWith(cstr{jc}, 'C')
                strt = cstr{jc};
                strt(1) = '!';
                cstr{jc} = strt;
            end
            if startsWith(cstr{jc}, '      ')
                strt = cstr{jc};
                strt = strt(7:end);
                cstr{jc} = strt;
            end
        end
        cstr(~startsWith(strtrim(cstr), '!') & ~contains(cstr, 'FORMAT')) = lower(cstr(~startsWith(strtrim(cstr), '!') & ~contains(cstr, 'FORMAT')));
        fid = fopen(filename, 'w');  % Open/create file for writing. Discard existing contents, if any.
        if fid == -1
            error('Cannot open file %s', filename);
        end
        fprintf(fid, '%s\n', cstr{:});
        fclose(fid);
    end
end
