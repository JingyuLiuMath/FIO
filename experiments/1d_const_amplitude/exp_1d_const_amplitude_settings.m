clear;
close all;

if ispc
    data_path = "./data/";
elseif isunix
    data_path = "/scratch/jyliu/FIO/1d_const_amplitude/data/";
end

originalPath = path;
addpath('../../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_1D(x, xi);
k_func = @(x, xi) exp_phi_func(x, xi);

if ispc
    p_list = (10 : 2 : 14)';
elseif isunix
    p_list = (10 : 2 : 18)';
end
num_n = length(p_list);

r_bf = 10;
tol_bf = 1e-8;

if ispc
    min_points = 64;
elseif isunix
    min_points = 256;
end

rank_func_tol = @(tol) 4 * log10(1 / tol);
tol_hss_list = [1e-3];
num_tol_hss = length(tol_hss_list);
tol_hss_display_list = ["10^{-3}"];

tol_cg = 1e-12;
if ispc
    maxit_cg = 50;
elseif isunix
    maxit_cg = 500;
end

num_sample = 256;

figure_prefix = "./figure/1d_const_amplitude";
