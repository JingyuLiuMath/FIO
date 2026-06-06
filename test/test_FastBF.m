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
tol_bf = 1e-6;

num_sample = 256;

%% K.
x = (0 : (N - 1))' / N;
xi = (-half_N : (half_N - 1))';

% K = exp_phi_func(x, xi);

%% BF.
fprintf("BF.\n");

tic;
[K_BF, ~] = fastBF(exp_phi_func, x, xi, r_bf, tol_bf);
t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

nnz_BF = my_nnz_bf(K_BF, "bf");
nnz_BF = nnz_BF.nnz_U + sum(nnz_BF.nnz_GTol) + nnz_BF.nnz_M + sum(nnz_BF.nnz_HTol) + nnz_BF.nnz_V;
nnz_dense = N^2;
ratio = nnz_BF / nnz_dense;
fprintf("  ratio: %.1e\n", ratio);

%%
f_ex = randn(N,1) + 1i * randn(N,1);
tic;
Kf = apply_fbf(K_BF, f_ex);
t_BF_apply = toc;
rel_err_BF = fbf_check(N, exp_phi_func, f_ex, x, xi, Kf, num_sample);
fprintf("  t_BF_apply: %.1e\n", t_BF_apply);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

% f_ex = randn(N,1) + 1i * randn(N,1);
% Kf = apply_fbf_adj(K_BF, f_ex);
% Kf_ex = K' * f_ex;
% rel_err_BF = norm(Kf_ex - Kf) / norm(Kf_ex);
% fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

%% Remove path.
path(originalPath);
