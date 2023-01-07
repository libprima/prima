function mod_files = list_mod_files(dir_name)
%LIST_MOD_FILES lists all module files (*.mod) in a directory 

mod_files = files_with_wildcard(dir_name, '*.mod');

return
