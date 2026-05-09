%% Setting.
clear;
close all;

originalPath = path;
addpath('../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_2D(x, xi);

p = 6;
n = 2^p;
N = n^2;
half_n = n / 2;

r_bf = 10;
tol_bf = 1e-8;

num_sample = 256;

min_points = 256;
tol_hss = 1e-3;

tol_cg = 1e-12;
maxit_cg = 500;

%% Run.
fprintf("\n");
result = run_FIO2D_inv_CG(...
    exp_phi_func, n, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg);

f_ex = result.f_ex;
Kf_ex = result.Kf_ex;
N = result.N;

fprintf("\n");
result = [];
result = run_FIO2D_inv(...
    exp_phi_func, n, ...
    r_bf, tol_bf, ...
    min_points, tol_hss, ...
    num_sample, ...
    tol_cg, maxit_cg, ...
    f_ex, Kf_ex);

path(originalPath);
