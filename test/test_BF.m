%% Setting.
clear;
close all;

addpath('../extern/BF.m/1D/src');
addpath('../extern/BF.m/1D/test');

c_func = @(x) (2 + sin(2 * pi * x)) / 8;
phi_func = @(x, xi) x * xi.' + c_func(x) * abs(xi).';

p = 10;
N = 2^p;
half_N = N / 2;

func_name = 'fun0';

r_BF = 8;
tol_BF = 1e-12;

num_sample = 256;

%% K.
x = (0 : (N - 1))' / N;
xi = (-half_N : (half_N - 1))';

K = exp(2 * pi * 1i * phi_func(x, xi));

%% BF.
fprintf("BF.\n");
BF_input = struct();
BF_input.xbox = [1,N+1];
BF_input.xx = (1:N)';
BF_input.kbox = [1,N+1];
BF_input.kk = (1:N)';
BF_input.fun = @(x,k)fun0(N,x,k);
BF_input.mR = r_BF;
BF_input.tol = tol_BF;

tic;
K_BF = bf_explicit(BF_input.fun, ...
    BF_input.xx, BF_input.xbox, ...
    BF_input.kk, BF_input.kbox, ...
    BF_input.mR, BF_input.tol, 0);
t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

f = randn(N,1) + 1i * randn(N,1);
Kf = apply_bf(K_BF, f);
rel_err_BF = bf_explicit_check(N, BF_input.fun, f, ...
    BF_input.xx, BF_input.kk, Kf, num_sample);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);
