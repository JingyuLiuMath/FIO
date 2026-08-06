clear;
close all;

if ispc
    data_path = "./data/";
elseif isunix
    data_path = ...
        "/scratch/jyliu/FIO/fio_inv_newton_schulz/2d_var_amplitude/data/";
end

if ~isfolder(data_path)
    mkdir(data_path);
end

a_func = @(x, xi) a_fun_2D_var(x, xi);
phi_func = @(x, xi) phi_fun_2D_var(x, xi);

p_list = (4 : 6)';
num_n = length(p_list);

n_leaf_bf = 8;
r_bf = 10;
tol_bf = 1e-8;
type_bf = "fbf";
type_bf_inv = "bf";

r_ns = 64;
tol_ns = 1e-1;
maxit_ns = 10;
num_test = 8;

restart_gmres = 50;
tol_gmres = 1e-12;
maxit_gmres = 10;

num_sample = 256;
