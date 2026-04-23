%% Setting.
clear;
close all;

addpath('../extern/FastBF.m/src');
addpath('../extern/FastBF.m/test/kernels');

exp_phi_func = @(x, xi) fun0_2D(x, xi);

p = 6;
n = 2^p;
N = n^2;
half_n = n / 2;

r_BF = 10;
tol_BF = 1e-6;

num_sample = 256;

%% K.
x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);
xi_co = (-half_n : (half_n - 1))';
xi = TensorProduct2D(xi_co, xi_co);

K = exp_phi_func(x, xi);

%% BF.
fprintf("BF.\n");

tic;
[K_BF, ~] = fastBF(exp_phi_func, x, xi, r_BF, tol_BF);
t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

f_ex = randn(N,1) + 1i * randn(N,1);
Kf = apply_fbf(K_BF, f_ex);
rel_err_BF = fbf_check(n, exp_phi_func, f_ex, x, xi, Kf, num_sample);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

f_ex = randn(N,1) + 1i * randn(N,1);
Kf = apply_fbf_adj(K_BF, f_ex);
Kf_ex = K' * f_ex;
rel_err_BF = norm(Kf_ex - Kf) / norm(Kf_ex);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);
