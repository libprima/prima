function iw = iswritable(fdname)
%ISWRITABLE returns true or false to indicate whether the file or directory is writable.
[status, attributes] = fileattrib(fdname);

iw = false;
if status == 1
    iw = attributes.UserWrite;
else
    errid = 'IsWritable:FDNotExists';
    error(errid, 'IsWritable:Error checking the attribute of %s. The file or directory does not exist', fdname);
end
