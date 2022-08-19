function hdir = homedir()
%This function gets the full path to the home directory of the current user.

hdir ='';
if isunix || ismac
    hdir = getenv('HOME');
elseif ispc
    hdir = getenv('HOMEPATH');
end

if isempty(hdir)
    try
        % According to
        % https://stackoverflow.com/questions/35887777/how-to-have-home-in-matlabs-save-saveas
        % The following should work on any platform with JVM.
        hdir = char(java.lang.System.getProperty('user.home'));
    catch
        error('HomeDir:Fail', 'Fail to get the path to the home directory.');
    end
end
