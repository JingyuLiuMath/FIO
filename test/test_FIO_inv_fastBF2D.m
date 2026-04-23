%% Setting.
clear;
close all;

addpath('../extern/FastBF.m/src');
addpath('../extern/FastBF.m/test/kernels');

exp_phi_func = @(x, xi) fun0_2D(x, xi);

p = 5;
n = 2^p;
N = n^2;
half_n = n / 2;

r_bf = 10;
tol_bf = 1e-8;

num_sample = 256;

min_points = 256;
r_hss = 2 * n;
tol_hss = 1e-3;

tol_cg = 1e-12;
maxit_cg = 500;

exact_flag = 0;

%% K
if exact_flag == 1
    x_co = (0 : (n - 1))' / n;
    x = TensorProduct2D(x_co, x_co);
    xi_co = (-half_n : (half_n - 1))';
    xi = TensorProduct2D(xi_co, xi_co);

    fprintf("exact info.\n")
    K = exp_phi_func(x, xi);
    r_K = rank(K);
    fprintf("  rank(K): %d\n", r_K);
    fprintf("  cond(K): %.1e\n", cond(K));
    if r_K ~= N
        keyboard;
    end

    xi_co1 = (-half_n : -1)';
    xi_co2 = (0 : (half_n - 1))';
    
    sub_xi_row = TensorProduct2D(xi_co1, xi_co1);
    sub_xi_col = TensorProduct2D(xi_co1, xi_co2);

    sub_K_row = exp_phi_func(x, sub_xi_row);
    sub_K_col = exp_phi_func(x, sub_xi_col);

    subG = sub_K_row' * sub_K_col;
    fprintf("  rank(subG, tol): %d\n", rank(subG, tol_hss * norm(subG)));

    % ShowSingularValueDecay(subG);
end

%% Run.
result = run_FIO_inv_fastBF2D(...
    exp_phi_func, n, ...
    r_bf, tol_bf, ...
    min_points, r_hss, tol_hss, ...
    num_sample, ...
    tol_cg, maxit_cg);
