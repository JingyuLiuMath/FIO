function fio_startup()
% fio_startup adds the project and required dependency paths.

root_path = fileparts(mfilename("fullpath"));

extern_path = fullfile(root_path, "extern");
if isfolder(extern_path)
    src_dir_list = dir(fullfile(extern_path, "**", "src"));
    src_dir_list = src_dir_list([src_dir_list.isdir]);
    src_path_list = strings(length(src_dir_list), 1);
    for i = 1 : length(src_dir_list)
        src_path_list(i) = fullfile(...
            src_dir_list(i).folder, src_dir_list(i).name);
    end
    src_path_list = sort(src_path_list);
    for i = 1 : length(src_path_list)
        addpath(genpath(src_path_list(i)), "-end");
    end
else
    warning("FIO:ExternNotFound", ...
        "The extern directory was not found at %s.", extern_path);
end

src_path = fullfile(root_path, "src");
addpath(genpath(src_path));

end
