function copy_shared_tools(setup_dir, runtime_dir)
%COPY_SHARED_TOOLS copies some tools from the setup directory to he runtime directory. These tools
% are shared between the setup and the runtime of the package. Their source code is maintained in
% the setup directory, and they need to be copied to the runtime directory during the setup.
%
% N.B.: Some of the tools are created during the setup, including `all_precisions.m` and
% `all_variants.m`. We should call the current function to copy them only after they are created!

tool_list = {'all_solvers.m', 'all_precisions.m', 'all_variants.m', 'dbgstr.m', 'get_mexname.m'};

cellfun(@(tool_name) copyfile(fullfile(setup_dir, tool_name), runtime_dir), tool_list);

return
