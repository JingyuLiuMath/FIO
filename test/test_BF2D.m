%% Setting.
clear;
close all;

a_func = @(x, xi) ones(size(x, 1), size(xi, 1));
phi_func = @(x, xi) phi_fun_2D(x, xi);
exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));
k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

p_lim = 6;

p = 6;
n = 2^p;
N = n^2;
half_n = n / 2;

n_leaf_bf = 0;
r_bf = 24;
tol_bf = 1e-6;
type_bf = "bf";

num_sample = 256;

%% K.
x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);
xi_co = (-half_n : (half_n - 1))';
xi = TensorProduct2D(xi_co, xi_co);

if p <= p_lim
    K = k_func(x, xi);
end

%% BF.
fprintf("BF.\n");

tic;
K_BF = my_construct_bf(...
    n, a_func, phi_func, ...
    x, xi, ...
    n_leaf_bf, r_bf, tol_bf, type_bf);
t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

nnz_BF = my_nnz_bf(K_BF, type_bf);
nnz_dense = N^2;
ratio = nnz_BF / nnz_dense;
fprintf("  ratio: %.1e\n", ratio);

%% Apply.
f_ex = randn(N, 1) + 1i * randn(N, 1);
tic;
Kf = my_apply_bf(K_BF, f_ex, type_bf);
t_BF_apply = toc;
rel_err_BF = my_check_bf(...
    N, k_func, f_ex, x, xi, Kf, num_sample);
fprintf("  t_BF_apply: %.1e\n", t_BF_apply);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

if p <= p_lim
    f_ex = randn(N, 1) + 1i * randn(N, 1);
    Kf = my_apply_bf_adj(K_BF, f_ex, type_bf);
    Kf_ex = K' * f_ex;
    rel_err_BF = norm(Kf_ex - Kf) / norm(Kf_ex);
    fprintf("  rel_err_BF: %.1e\n", rel_err_BF);
end
