clear;
close all;

if ispc
    data_path = "./data/";
elseif isunix
    data_path = "/scratch/jyliu/FIO/2d_const_amplitude/data/";
end

originalPath = path;
addpath('../../extern/FastBF.m/src');

a_func = @(x, xi) ones(size(x, 1), size(xi, 1));
phi_func = @(x, xi) phi_fun_2D(x, xi);
exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));
k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

if ispc
    p_list = (4 : 6)';
elseif isunix
    p_list = (6 : 9)';
end
num_n = length(p_list);

n_leaf_bf = 8;
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
rank_func_tol = @(tol) ceil(3 * log10(1 / tol));

tol_cg = 1e-12;
if ispc
    maxit_cg = 50;
elseif isunix
    maxit_cg = 500;
end

num_sample = 256;

figure_prefix = "./figure/2d_const_amplitude";
