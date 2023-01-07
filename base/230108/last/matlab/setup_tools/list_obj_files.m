function obj_files = list_obj_files(dir_name)
%LIST_OBJ_FILES lists all object files (*.o, *.obj) in a directory

obj_files = [files_with_wildcard(dir_name, '*.o'), files_with_wildcard(dir_name, '*.obj')];

return
