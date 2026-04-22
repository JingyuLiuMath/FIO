%% Setting.
clear;
close all;

addpath('../extern/FastBF.m/src');

c_func = @(x) (2 + sin(2 * pi * x)) / 8;
phi_func = @(x, xi) x * xi.' + c_func(x) * abs(xi).';
exp_phi_func = @(x, xi) exp(2 * pi * 1i * phi_func(x, xi));

p = 10;
N = 2^p;
half_N = N / 2;

r_bf = 10;
tol_bf = 1e-8;

num_sample = 256;

min_points = 256;
r_hss = 2 * p;
tol_hss = 1e-3;

tol_cg = 1e-12;
maxit_cg = 500;

exact_flag = 1;

%% K
if exact_flag == 1
    x = (0 : (N - 1))' / N;
    xi = (-half_N : (half_N - 1))';

    fprintf("exact info.\n")
    K = exp_phi_func(x, xi);
    r_K = rank(K);
    fprintf("  rank(K): %d\n", r_K);
    fprintf("  cond(K): %.1e\n", cond(K));
    if r_K ~= N
        keyboard;
    end

    G = K' * K;

    subG = G(1 : half_N, (half_N + 1) : end);
    fprintf("  rank(subG, tol): %d\n", rank(subG, tol_hss * norm(subG)));

    % ShowSingularValueDecay(subG);
end

%% Run.
result = run_FIO_inv_fastBF(...
    exp_phi_func, N, ...
    r_bf, tol_bf, ...
    min_points, r_hss, tol_hss, ...
    num_sample, ...
    tol_cg, maxit_cg);
