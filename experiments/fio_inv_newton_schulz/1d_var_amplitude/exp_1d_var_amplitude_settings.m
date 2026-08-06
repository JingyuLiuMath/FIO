clear;
close all;

if ispc
    data_path = "./data/";
elseif isunix
    data_path = ...
        "/scratch/jyliu/FIO/fio_inv_newton_schulz/1d_var_amplitude/data/";
end

if ~isfolder(data_path)
    mkdir(data_path);
end

originalPath = path;

phi_func = @(x, xi) phi_fun_1D(x, xi);

m = 10;
sigma_sq = 0.1;
rng(0);
x_pts = rand(m, 1);
xi_pts = rand(m, 1) - 0.5;

p_list = (10 : 2 : 14)';
num_n = length(p_list);

if ispc
    n_leaf_bf = 8;
else
    n_leaf_bf = 16;
end
r_bf = 10;
tol_bf = 1e-8;
type_bf = "mybf";
type_bf_inv = "mybf";
type_bf_smp = "fbf";
type_bf_inv_smp = "bf";

r_ns = 4;
tol_ns = 4e-1;
maxit_ns = 20;
num_test = 8;

restart_gmres = 50;
tol_gmres = 1e-12;
if ispc
    maxit_gmres = 10;
elseif isunix
    maxit_gmres = 50;
end

num_sample = 256;
