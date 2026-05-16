%% Setting.
clear;
close all;

originalPath = path;
addpath('../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_2D(x, xi);

p = 5;
n = 2^p;
N = n^2;

r_bf = 10;
tol_bf = 1e-8;

num_sample = 256;

min_points = 64;
tol_hss = 1e-3;

tol_cg = 1e-12;
maxit_cg = 50;

%% HSS Rank.
% fprintf("\n");
% result_rank = run_FIO2D_hss_rank(...
%     exp_phi_func, n);

%% BF.
fprintf("\n");
result_bf = run_FIO2D_inv_CG(...
    exp_phi_func, n, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg);

%% FIO.
fprintf("\n");
run_FIO2D_inv(...
    result_bf, ...
    min_points, tol_hss);

path(originalPath);
