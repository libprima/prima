


function str = str_fixed(str, first, last, intent)

spaces = '                                                                ';

str = regexprep(str, '\t', '    ');

first_non_space = min(find(~isspace(str), 1, 'first'));
num_leading_spaces = first - 1 + first_non_space - 1;

leading_spaces = spaces(1 : num_leading_spaces);

col = last - num_leading_spaces;
str = strtrim(str);
len = length(str);

if (len <= col)
    row = 1;
else
    row = ceil((len - col)/(col - 2)) + 1;
end

strnew = [leading_spaces, str(1 : min(len, col))];

for i = 2 : row
    str_tmp = [leading_spaces, spaces(1:intent), str(col + (i-2)*(col - 2) + 1 : min(len, col + (i-1)*(col - 2)))];
    str_tmp(first - 1) = '&';
    str_tmp = ['\n', str_tmp];
    strnew = [strnew, str_tmp];
end

str = strnew;
