function ndftpath = showpath()

fprintf('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n');
fprintf('\nCurrently in %s\n', pwd());

fprintf('\nNon-default paths:\n');
oldpath = path();
restoredefaultpath;
dftpath = path();
path(oldpath);
ndftpath = setdiff(split(oldpath, ':'), split(dftpath, ':'));
for ic = 1 : length(ndftpath)
    fprintf('\n%s', ndftpath{ic});
end

fprintf('\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n');
