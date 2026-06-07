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

p_lim = 12;

p = 10;
N = 2^p;
half_N = N / 2;

r_bf = 10;
tol_bf = 1e-6;

num_sample = 256;

%% K.
x = (0 : (N - 1))' / N;
xi = (-half_N : (half_N - 1))';

if p <= p_lim
    K = k_func(x, xi);
end

%% BF.
fprintf("BF.\n");

tic;
[K_BF, ~] = fastBF(k_func, x, xi, r_bf, tol_bf);
t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

nnz_BF = my_nnz_bf(K_BF, "bf");
nnz_BF = nnz_BF.nnz_U + sum(nnz_BF.nnz_GTol) + nnz_BF.nnz_M + sum(nnz_BF.nnz_HTol) + nnz_BF.nnz_V;
nnz_dense = N^2;
ratio = nnz_BF / nnz_dense;
fprintf("  ratio: %.1e\n", ratio);

%% Apply.
f_ex = randn(N, 1) + 1i * randn(N, 1);
tic;
Kf = apply_fbf(K_BF, f_ex);
t_BF_apply = toc;
rel_err_BF = fbf_check(N, k_func, f_ex, x, xi, Kf, num_sample);
fprintf("  t_BF_apply: %.1e\n", t_BF_apply);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

if p <= p_lim
    f_ex = randn(N, 1) + 1i * randn(N, 1);
    Kf = apply_fbf_adj(K_BF, f_ex);
    Kf_ex = K' * f_ex;
    rel_err_BF = norm(Kf_ex - Kf) / norm(Kf_ex);
    fprintf("  rel_err_BF: %.1e\n", rel_err_BF);
end

%% Remove path.
path(originalPath);
