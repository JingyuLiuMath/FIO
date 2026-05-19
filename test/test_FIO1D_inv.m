%% Setting.
clear;
close all;

originalPath = path;
addpath('../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_1D(x, xi);

p = 10;
N = 2^p;
half_N = N / 2;

r_bf = 10;
tol_bf = 1e-8;

num_sample = 256;

min_points = 64;
r_hss = 25;
tol_hss = 1e-3;

tol_cg = 1e-12;
maxit_cg = 50;

indep = 0;

%% HSS Rank.
fprintf("\n");
result_rank = run_FIO1D_hss_rank(...
    exp_phi_func, N);

%% BF.
fprintf("\n");
result_bf = run_FIO1D_inv_CG(...
    exp_phi_func, N, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg);

%% FIO.
fprintf("\n");
run_FIO1D_inv(...
    result_bf, ...
    min_points, tol_hss, indep);

path(originalPath);
