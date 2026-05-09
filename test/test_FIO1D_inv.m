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

min_points = 256;
r_hss = 25;
tol_hss = 1e-3;

tol_cg = 1e-12;
maxit_cg = 500;

%% Run.
fprintf("\n");
result = run_FIO1D_inv_CG(...
    exp_phi_func, N, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg);

f_ex = result.f_ex;
Kf_ex = result.Kf_ex;
N = result.N;

fprintf("\n");
result = [];
result = run_FIO1D_inv(...
    exp_phi_func, N, ...
    r_bf, tol_bf, ...
    min_points, tol_hss, ...
    num_sample, ...
    tol_cg, maxit_cg, ...
    f_ex, Kf_ex);

path(originalPath);
