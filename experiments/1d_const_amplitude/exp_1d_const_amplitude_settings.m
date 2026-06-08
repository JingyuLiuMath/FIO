clear;
close all;

if ispc
    data_path = "./data/";
elseif isunix
    data_path = "/scratch/jyliu/FIO/1d_const_amplitude/data/";
end

originalPath = path;
addpath('../../extern/FastBF.m/src');

a_func = @(x, xi) ones(size(x, 1), size(xi, 1));
phi_func = @(x, xi) phi_fun_1D(x, xi);
exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));
k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

if ispc
    p_list = (10 : 2 : 14)';
elseif isunix
    p_list = (10 : 2 : 18)';
end
num_n = length(p_list);

if ispc
    n_leaf_bf = 8;
else
    n_leaf_bf = 16;
end
r_bf = 10;
tol_bf = 1e-8;
type_bf = "bf";

if ispc
    N_leaf_hss = 64;
elseif isunix
    N_leaf_hss = 256;
end
tol_hss_list = [1e-3];
num_tol_hss = length(tol_hss_list);
tol_hss_display_list = ["10^{-3}"];
rank_func_tol = @(tol) 4 * log10(1 / tol);

tol_cg = 1e-12;
if ispc
    maxit_cg = 50;
elseif isunix
    maxit_cg = 500;
end

num_sample = 256;

figure_prefix = "./figure/1d_const_amplitude";
