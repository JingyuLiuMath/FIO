function fio_startup()
% fio_startup adds the project and required dependency paths.

root_path = fileparts(mfilename("fullpath"));

fastbf_path = fullfile(root_path, "extern", "FastBF.m", "src");
if isfolder(fastbf_path)
    addpath(fastbf_path);
else
    warning("FIO:FastBFNotFound", ...
        "FastBF was not found at %s.", fastbf_path);
end

src_path = fullfile(root_path, "src");
addpath(genpath(src_path));

end