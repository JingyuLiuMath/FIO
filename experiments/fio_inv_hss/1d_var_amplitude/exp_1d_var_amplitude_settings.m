clear;
close all;

if ispc
    data_path = "./data/";
elseif isunix
    data_path = "/scratch/jyliu/FIO/fio_inv_hss/1d_var_amplitude/data/";
end

originalPath = path;

phi_func = @(x, xi) phi_fun_1D(x, xi);
exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));

m = 10;
sigma_sq = 0.1;
x_pts = rand(m, 1);
xi_pts = (rand(m, 1) - 0.5);

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
type_bf = "fbf";

if ispc
    N_leaf_hss = 64;
elseif isunix
    N_leaf_hss = 256;
end

tol_hss_list = [1e-3];
num_tol_hss = length(tol_hss_list);
tol_hss_display_list = ["10^{-3}"];
rank_func_tol = @(tol) 8 * log10(1 / tol);

tol_cg = 1e-12;
if ispc
    maxit_cg = 50;
elseif isunix
    maxit_cg = 500;
end

num_sample = 256;

figure_prefix = "./figure/1d_var_amplitude";
