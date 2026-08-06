%% Setting.
clear;
close all;

a_func = @(x, xi) ones(size(x, 1), size(xi, 1));
phi_func = @(x, xi) phi_fun_2D(x, xi);
exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));
k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

p = 5;
n = 2^p;
N = n^2;

n_leaf_bf = 8;
r_bf = 10;
tol_bf = 1e-8;
type_bf = "fbf";
type_bf_inv = "bf";

num_sample = 256;

restart_gmres = 50;
tol_gmres = 1e-12;
maxit_gmres = 10;

r_ns = 8;
tol_ns = 1e-1;
maxit_ns = 10;
num_test = 8;

%% BF.
fprintf("\n");
result_bf = run_FIO2D_inv_GMRES(...
    n, ...
    a_func, phi_func, ...
    n_leaf_bf, r_bf, tol_bf, type_bf, ...
    num_sample, ...
    restart_gmres, tol_gmres, maxit_gmres);

%% Newton-Schulz.
fprintf("\n");
run_FIO2D_inv_NewtonSchulz(...
    result_bf, ...
    r_ns, tol_ns, maxit_ns, num_test, type_bf_inv);
