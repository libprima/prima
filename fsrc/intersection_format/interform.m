function interform(directory)
%INTERFORM refactors the Fortran files in directory into the "intersection
% format" and save them in outputdir.
% See http://fortranwiki.org/fortran/show/Continuation+lines for details.
% 
% Coded by Zaikun Zhang in August, 2020.

outputdir = 'intersection_format';

if nargin < 1
    directory = cd();
end
cd(directory);
givendir = cd();  % Full path of the given directory, which is the current directory now.

outputdir = fullfile(givendir, outputdir);  % Full path of the output directory.
if exist(outputdir, 'dir')
    rmdir(outputdir, 's');
end
mkdir(outputdir);  % Make the output directory
copyfile([mfilename('fullpath'), '.m'], outputdir);  % Save the current script in the output directory

readme = 'README.txt';
fid = fopen(fullfile(outputdir, readme), 'w');
if fid == -1
    error('Cannot open file %s.', readme);
end
fprintf(fid, 'This folder contains the intersection-format version of the Fortran source\n');
fprintf(fid, 'files. The files in this folder are generated automatically by the %s.m\n', mfilename); 
fprintf(fid, 'script, and they are NOT intended to be readable.\n');
fprintf(fid, '\n');
fprintf(fid, 'This project is coded in the free format, yet some platforms accept only\n');
fprintf(fid, 'fixed-format Fortran code, for example, the MATLAB MEX on Windows. The code\n');
fprintf(fid, 'in this folder can serve such a purpose.\n');
fprintf(fid, '\n');
fprintf(fid, 'In the intersection format, each continued line has an ampersand at\n'); 
fprintf(fid, 'column 73, and each continuation line has an ampersand at column 6.\n');
fprintf(fid, 'A Fortran file in such a format can be compiled both as a fixed-format\n');
fprintf(fid, 'file and as a free-format file.\n');
fprintf(fid, 'See http://fortranwiki.org/fortran/show/Continuation+lines for details.\n');
fprintf(fid, '\n');
fprintf(fid, 'Zaikun Zhang (www.zhangzk.net), %s', date);
fclose(fid);

% Do not perform refactoring in these subdirectories (if exist)
ignoredir = {'original', 'classical', 'backup', 'intersection_format', 'trash', 'test'};  

% Ignore the following files
ignorefile = {'calfun__genmod.f90', 'mexfunction__genmod.f90', 'test.f'};


% The following lines perform the refactoring in the current directory
ffiles = [dir(fullfile(givendir, '*.f90')); dir(fullfile(givendir, '*.F90')); dir(fullfile(givendir, '*.f')); dir(fullfile(givendir, '*.F'))];
ffiles = setdiff({ffiles.name}, ignorefile);
for j = 1 : length(ffiles)
    copyfile(fullfile(givendir, ffiles{j}), outputdir);   
    [~, ~, ext] = fileparts(ffiles{j});
    if strcmpi(ext, '.f90')
        refactor_file(fullfile(outputdir, ffiles{j}));
    end
end
delete(fullfile(outputdir, '*.bak'));
delete(fullfile(outputdir, '*.f90'));
delete(fullfile(outputdir, '*.F90'));

% Copy the filelist to the output directory.
if exist(fullfile(givendir, 'filelist'), 'file')
    copyfile(fullfile(givendir, 'filelist'), outputdir);
    refactor_filelist(fullfile(outputdir, 'filelist'));
end

% Copy the header files to the output directory.
hfiles = dir(fullfile(givendir, '*.h')); 
hfiles = {hfiles.name};
for j = 1 : length(hfiles)
    copyfile(fullfile(givendir, hfiles{j}), outputdir);
end

% The following lines get a cell array containing the names (but not
% full path) of all the subdirectories of the given directory.
d = dir(givendir);
isub = [d(:).isdir]; 
subdir = {d(isub).name};
subdir = setdiff(subdir, [{'.','..'}, ignoredir]);

% The following lines perform the refactoring in the subdirectories of
% the current directory.
for i = 1 : length(subdir)
    ffiles = [dir(fullfile(givendir, subdir{i},'*.f90')); dir(fullfile(givendir, subdir{i}, '*.F90')); dir(fullfile(givendir, subdir{i}, '*.f')); dir(fullfile(givendir, subdir{i}, '*.F'))];
    ffiles = setdiff({ffiles.name}, ignorefile);
    
    mkdir(fullfile(outputdir, subdir{i}));
    for j = 1 : length(ffiles)
        copyfile(fullfile(givendir, subdir{i}, ffiles{j}), fullfile(outputdir, subdir{i}));   
        [~, ~, ext] = fileparts(ffiles{j});
        if strcmpi(ext, '.f90')
            refactor_file(fullfile(outputdir, subdir{i}, ffiles{j}));
        end
    end
    delete(fullfile(outputdir, subdir{i}, '*.bak'));
    delete(fullfile(outputdir, subdir{i}, '*.f90'));
    delete(fullfile(outputdir, subdir{i}, '*.F90'));
              
    % Copy the header files to the output directory.
    hfiles = dir(fullfile(givendir, subdir{i}, '*.h')); 
    hfiles = {hfiles.name};
    for j = 1 : length(hfiles)
        copyfile(fullfile(givendir, subdir{i}, hfiles{j}), fullfile(outputdir, subdir{i}));
    end

    % Copy the filelist to the output directory.
    if exist(fullfile(givendir, subdir{i}, 'filelist'), 'file')
        copyfile(fullfile(givendir, subdir{i}, 'filelist'), fullfile(outputdir, subdir{i}));
        refactor_filelist(fullfile(outputdir, subdir{i}, 'filelist'));
    end
    
end


function refactor_file(filename)
%REFACTOR_FILE refactors a given file into the "intersection format". 
% See http://fortranwiki.org/fortran/show/Continuation+lines for details.
% The new file has the same name as the original file, but the extension
% ".f90" will be changed to ".f", and ".F90" will be changed to ".F". If
% the new file has the same name as the original one, then the original
% file will be backuped in "ORIGINAL_FILE_NAME.bak".

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file %s.', filename);
end

% Read the file into a cell of strings 
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);
cstr = data{1};
cstr = deblank(cstr);  % Remove trailing blanks

i = 1;
j = 1;
k = 1;

while(j <= length(cstr))
    strtmp = cstr{j};
    strtmp_trimed = strtrim(strtmp);
    
    if isempty(strtmp_trimed) || strcmp(strtmp_trimed(1), '!') || strcmp(strtmp_trimed(1), '#')
        strs{k} = strtmp_trimed;
        k = k + 1;
        i = j + 1;
        j = j + 1;
    elseif j < length(cstr) && ~isempty(strtmp) && strcmp(strtmp(end), '&') 
        cstr{j} = strtmp(1 : end-1);
        strtmp = strtrim(cstr{j+1});
        if (~isempty(strtmp) && strcmp(strtmp(1), '&'))
            strtmp = strtmp(2 : end);
            cstr{j+1} = strtmp;
        end
        j = j + 1;
    else
        strs{k} = strjoin(cstr(i:j), ' ');
        strs{k} = refactor_str(strs{k}, 7, 72);
        k = k + 1;
        i = j + 1;
        j = j + 1;
    end
end 

% Save the refactored file 
refactored_filename = regexprep(filename, '.f90$', '.f');
refactored_filename = regexprep(refactored_filename, '.F90$', '.F');
if strcmp(refactored_filename, filename)
    copyfile(filename, [filename, '.bak']);
end
fid = fopen(refactored_filename, 'wt');
if fid == -1
    error('Cannot open file %s.', refactored_filename);
end

[~, fname, ext] = fileparts(filename);
fprintf(fid, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
fprintf(fid, '! This is the intersection-format version of %s.\n', [fname, ext]);
fprintf(fid, '! The file is generated automatically and is NOT intended to be readable.\n');
fprintf(fid, '!\n');
fprintf(fid, '! In the intersection format, each continued line has an ampersand at\n'); 
fprintf(fid, '! column 73, and each continuation line has an ampersand at column 6.\n');
fprintf(fid, '! A Fortran file in such a format can be compiled both as a fixed-format\n');
fprintf(fid, '! file and as a free-format file.\n');
fprintf(fid, '! See http://fortranwiki.org/fortran/show/Continuation+lines for details.\n');
fprintf(fid, '!\n');
fprintf(fid, '! Generated using the %s.m script by Zaikun Zhang (www.zhangzk.net)\n! on %s.\n', mfilename, date);
fprintf(fid, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n');

for i = 1 : length(strs)
    fprintf(fid, strs{i});
    if i < length(strs)
        fprintf(fid, '\n');
    end
end
fclose(fid);


function str = refactor_str(str, first, last)
%REFACTOR_STR refactors a given string into the "intersection format".
% See http://fortranwiki.org/fortran/show/Continuation+lines for details.

spaces = '                                                                ';

str = regexprep(str, '\t', '    ');  % Replace each tab by four spaces

first_non_space = min(find(~isspace(str), 1, 'first'));  % Index of the first non-space charactor in str
num_leading_spaces = first - 1 + first_non_space - 1;  % Number of leading spaces in str after refactoring

leading_spaces = spaces(1 : num_leading_spaces); % Save the leading spaces in a string

width_first_row = last - num_leading_spaces; % Width of the first row of str after refactoring
width = last - first + 1;  % Width of the other rows of str after refactoring

str = strtrim(str);  % Remove the leading and trailing spaces from str
str = regexprep(str,' +',' ');  % Replace all the continuous multiple spaces by one single space
len = length(str);  % Length of the trimed str

row = ceil((len - width_first_row)/width) + 1;  % Number of rows of str after refactoring

strnew = [leading_spaces, str(1 : min(len, width_first_row))];  % The first row after refactoring

for i = 2 : row
    strnew = [strnew, '&'];  % Append an '&' at the end of the i-1 th row.
    strtmp = strtrim(str(width_first_row + (i-2)*width + 1 : min(len, width_first_row + (i-1)*width)));  % Content of the i th row
    strtmp = [spaces(1:first-2), '&', strtmp];  % Add first - 2 spaces and an '&' at the beginning of the i-th row
    strtmp = ['\n', strtmp];  % Add a '\n' at the beginning of the i-th row
    strnew = [strnew, strtmp];
end

str = strnew;


function refactor_filelist(filename)

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file %s.', filename);
end

% Read the file into a cell of strings 
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);
cstr = data{1};
cstr = deblank(cstr);  % Remove trailing blanks

for i = 1 : length(cstr)
    cstr{i} = regexprep(cstr{i}, '.F90$', '.F');
    cstr{i} = regexprep(cstr{i}, '.f90$', '.f');
end

% Save the file again
fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file %s.', filename);
end
fprintf(fid, '%s\n', cstr{:});
fclose(fid);
