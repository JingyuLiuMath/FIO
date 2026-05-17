clear;
close all;

profile('-memory','on');  % an alternative

if ispc
    data_path = "./data/";
else
    data_path = "/scratch/jyliu/FIO/2d_const_amplitude/data/";
end

originalPath = path;
addpath('../../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_2D(x, xi);

if ispc
    p_list = (4 : 6)';
else
    p_list = (6 : 9)';
end
num_n = length(p_list);

r_bf = 10;
tol_bf = 1e-8;

if ispc
    min_points = 64;
else
    min_points = 256;
end
tol_hss_list = [1e-3];
num_tol_hss = length(tol_hss_list);

tol_cg = 1e-12;
if ispc
    maxit_cg = 50;
else
    maxit_cg = 500;
end

num_sample = 256;