%% Setting.
clear;
close all;

addpath('../extern/BF.m/2D/GBF/src');
addpath('../extern/BF.m/2D/MBF/src');

exp_phi_func = @(x, xi) fun_2D(x, xi);

p = 5;
n = 2^p;
N = n^2;
half_n = n / 2;

r_bf = 10;
tol_bf = 1e-6;

num_sample = 256;

%% K.
x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);
xi_co = (-half_n : (half_n - 1))';
xi = TensorProduct2D(xi_co, xi_co);

K = exp_phi_func(x, xi);

xbox = [0, 1; 0, 1];
kbox = [-half_n, half_n; -half_n, half_n];

%% BF.
fprintf("BF.\n");

tic;
K_BF = mbf_explicit(exp_phi_func, x, xbox, xi, kbox, r_bf, tol_bf, 0);
t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

f_ex = randn(N,1) + 1i * randn(N,1);
Kf = apply_mbf(K_BF, f_ex);
rel_err_BF = mbf_explicit_check(n, exp_phi_func, f_ex, Kf, num_sample);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

f_ex = randn(N,1) + 1i * randn(N,1);
Kf = apply_mbf_adj(K_BF, f_ex);
Kf_ex = K' * f_ex;
rel_err_BF = norm(Kf_ex - Kf) / norm(Kf_ex);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

%% Remove path.
rmpath('../extern/BF.m/2D/GBF/src');
rmpath('../extern/BF.m/2D/MBF/src');
