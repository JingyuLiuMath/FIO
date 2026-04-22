%% Setting.
clear;
close all;

c_func = @(x) (2 + sin(2 * pi * x)) / 8;
phi_func = @(x, xi) x * xi.' + c_func(x) * abs(xi).';

p = 10;
N = 2^p;
half_N = N / 2;

func_name = 'fun0';

r_bf = 8;
tol_bf = 1e-12;

num_sample = 256;

min_points = 256;
r_hss = 2 * p;
tol_hss = 1e-3;

tol_cg = 1e-12;
maxit_cg = 500;

exact_flag = 1;

%% K
if exact_flag == 1
    fprintf("exact info.\n")
    x = (0 : (N - 1))' / N;
    xi = (-half_N : (half_N - 1))';

    K = exp(2 * pi * 1i * phi_func(x, xi));
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

%% BF
fprintf("BF.\n");
BF_input = struct();
BF_input.xbox = [1,N+1];
BF_input.xx = (1:N)';
BF_input.kbox = [1,N+1];
BF_input.kk = (1:N)';
BF_input.fun = @(x,k)fun0(N,x,k);
BF_input.mR = r_bf;
BF_input.tol = tol_bf;

tic;
K_BF = bf_explicit(BF_input.fun, ...
    BF_input.xx, BF_input.xbox, ...
    BF_input.kk, BF_input.kbox, ...
    BF_input.mR, BF_input.tol, 0);
t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

f_ex = randn(N,1) + 1i*randn(N,1);
tic;
Kf = apply_bf(K_BF, f_ex);
t_BF_apply = toc;

fprintf("  t_BF_apply: %.1e\n", t_BF_apply);
rel_err_BF = bf_explicit_check(N, BF_input.fun, f_ex, ...
    BF_input.xx, BF_input.kk, Kf, num_sample);
fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

%% HSS.
fprintf("HSS.\n");
op_G = @(v) apply_bf_adj(K_BF, apply_bf(K_BF, v));

tic;
G_HSS = BF_HSS(N);
G_HSS.BuildTree(min_points);
leaf_size = G_HSS.MaxLeafSize();
G_HSS.BlackBoxConstruct(op_G, leaf_size, r_hss, tol_hss);
t_HSS_construct = toc;
fprintf("  t_HSS_construct: %.1e\n", t_HSS_construct);

hss_rank = G_HSS.Rank();
fprintf("  HSS rank: %d\n", hss_rank);
hss_mem = G_HSS.Storage();
ratio = hss_mem / N^2;
fprintf("  ratio: %.1e\n", ratio);

tic;
G_HSS.ULV_Factor();
t_HSS_factor = toc;
fprintf("  t_HSS_factor: %.1e\n", t_HSS_factor);

%% Direct Solution.
fprintf("Direct solution.\n");
tic;
f_direct = G_HSS.ULV_Solve(apply_bf_adj(K_BF, Kf));
t_solve_direct = toc;

fprintf("  t_solve_direct: %.1e\n", t_solve_direct);
rel_res_direct = norm(Kf - apply_bf(K_BF, f_direct)) / norm(Kf);
fprintf("  rel_res_direct: %.1e\n", rel_res_direct);
rel_err_direct = norm(f_ex - f_direct) / norm(f_ex);
fprintf("  rel_res_direct: %.1e\n", rel_err_direct);

%% Iterative solution.
fprintf("Iterative solution.\n");
rhs = apply_bf_adj(K_BF, Kf);

fprintf("  Without precond.\n");
tic;
[f_cg, flag_cg, ~, iter_cg] = pcg(op_G, rhs, tol_cg, maxit_cg);
t_cg = toc;

fprintf("    t_cg: %.1e\n", t_cg);
fprintf("    iter_cg: %d\n", iter_cg);
rel_res_cg = norm(Kf - apply_bf(K_BF, f_cg)) / norm(Kf);
fprintf("    rel_res_cg: %.1e\n", rel_res_cg);
rel_err_cg = norm(f_ex - f_cg) / norm(f_ex);
fprintf("    rel_err_cg: %.1e\n", rel_err_cg);

fprintf("  With precond.\n");
op_M = @(v) G_HSS.ULV_Solve(v);
tic;
[f_pcg, flag_pcg, ~, iter_pcg] = pcg(op_G, rhs, tol_cg, maxit_cg, op_M);
t_pcg = toc;

fprintf("    t_pcg: %.1e\n", t_pcg);
fprintf("    iter_pcg:%d\n", iter_pcg);
rel_res_pcg = norm(Kf - apply_bf(K_BF, f_pcg)) / norm(Kf);
fprintf("    rel_res_pcg: %.1e\n", rel_res_pcg);
rel_err_pcg = norm(f_ex - f_pcg) / norm(f_ex);
fprintf("    rel_err_pcg: %.1e\n", rel_err_pcg);
