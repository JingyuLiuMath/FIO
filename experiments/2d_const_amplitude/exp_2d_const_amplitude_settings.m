clear;
close all;

if ispc
    data_path = "./data/";
elseif isunix
    data_path = "/scratch/jyliu/FIO/2d_const_amplitude/data/";
end

originalPath = path;
addpath('../../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_2D(x, xi);

if ispc
    p_list = (4 : 6)';
elseif isunix
    p_list = (6 : 9)';
end
num_n = length(p_list);

r_bf = 10;
tol_bf = 1e-8;

if ispc
    min_points = 64;
elseif isunix
    min_points = 256;
end
tol_hss_list = [1e-3];
num_tol_hss = length(tol_hss_list);

tol_cg = 1e-12;
if ispc
    maxit_cg = 50;
elseif isunix
    maxit_cg = 500;
end

num_sample = 256;