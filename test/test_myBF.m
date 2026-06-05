%% Setting.
clear;
close all;

originalPath = path;
addpath('../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_1D(x, xi);
phi_func = @(x, xi) phi_fun_1D(x, xi);
a_func = @(x, xi) ones(size(x, 1), size(xi, 1));

p = 10;
N = 2^p;
half_N = N / 2;

r_bf = 10;
tol_bf = 1e-6;

num_sample = 256;

%% K.
x = (0 : (N - 1))' / N;
xi = (-half_N : (half_N - 1))';

%% BF.
fprintf("BF.\n");

tic;
K_BF = Butterfly(N, a_func, phi_func, r_bf, tol_bf);
t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

nnz_BF = K_BF.Nnz();
nnz_dense = N^2;
ratio = nnz_BF / nnz_dense;
fprintf("  ratio: %.1e\n", ratio);

f_ex = randn(N,1) + 1i * randn(N,1);
Kf = K_BF.Apply(f_ex);
rel_err_BF = fbf_check(N, exp_phi_func, f_ex, x, xi, Kf, num_sample);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);


% f_ex = randn(N,1) + 1i * randn(N,1);
% Kf =  K_BF.ApplyAdj(f_ex);
% Kf_ex = K' * f_ex;
% rel_err_BF = norm(Kf_ex - Kf) / norm(Kf_ex);
% fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

%% Remove path.
path(originalPath);
