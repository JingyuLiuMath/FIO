%% Setting.
clear;
close all;

originalPath = path;
addpath('../extern/FastBF.m/src');

a_func = @(x, xi) ones(size(x, 1), size(xi, 1));
phi_func = @(x, xi) phi_fun_1D(x, xi);
exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));
k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

p = 10;
N = 2^p;

n_leaf_bf = 8;
r_bf = 10;
tol_bf = 1e-8;
type_bf = "mybf";

num_sample = 256;

N_leaf_hss = 64;
tol_hss = 1e-3;
rank_func_tol = @(tol) 4 * log10(1 / tol);

tol_cg = 1e-12;
maxit_cg = 50;

%% HSS Rank.
fprintf("\n");
result_rank = run_FIO1D_hss_rank(...
    k_func, N, rank_func_tol);

%% BF.
fprintf("\n");
result_bf = run_FIO1D_inv_CG(...
    N, ...
    a_func, phi_func, ...
    n_leaf_bf, r_bf, tol_bf, type_bf, ...
    num_sample, ...
    tol_cg, maxit_cg);

%% FIO.
fprintf("\n");
run_FIO1D_inv_HSS(...
    result_bf, ...
    N_leaf_hss, rank_func_tol, tol_hss);

path(originalPath);
