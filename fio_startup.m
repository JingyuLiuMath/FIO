function fio_startup()
% fio_startup

file_path = mfilename('fullpath');
tmp = strfind(file_path, 'fio');
file_path = file_path(1:(tmp(end)-1));
addpath(genpath([file_path 'src']));
addpath(genpath([file_path 'extern']));

end