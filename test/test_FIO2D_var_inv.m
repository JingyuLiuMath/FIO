%% Setting.
clear;
close all;

originalPath = path;
addpath('../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_2D_var(x, xi);
a_func = @(x, xi) a_fun_2D_var(x, xi);

p = 5;
n = 2^p;
N = n^2;

r_a = 20;
tol_a = 1e-7;

r_bf = 10;
tol_bf = 1e-8;

num_sample = 256;

min_points = 64;
tol_hss = 1e-3;

tol_cg = 1e-12;
maxit_cg = 50;

indep = 1;


%% HSS Rank.
fprintf("\n");
result_rank = run_FIO2D_var_hss_rank(...
    exp_phi_func, n, a_func);

%% BF.
fprintf("\n");
result_bf = run_FIO2D_var_inv_CG(...
    a_func, r_a, tol_a, ...
    exp_phi_func, n, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg);

%% FIO.
fprintf("\n");
run_FIO2D_var_inv(...
    result_bf, ...
    min_points, tol_hss, indep);

path(originalPath);
