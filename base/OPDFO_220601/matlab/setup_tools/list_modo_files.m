function modo_files = list_modo_files(dir_name)
%LIST_MODO_FILES lists all module or object files in a directory

modo_files = [list_mod_files(dir_name), list_obj_files(dir_name)];

return
