function result = run_FIO_inv_fastBF2D(...
    exp_phi_func, n, ...
    r_bf, tol_bf, ...
    min_points, r_hss, tol_hss, ...
    num_sample, ...
    tol_cg, maxit_cg)

N = n^2;
fprintf("Basic info.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);

fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);

fprintf("  min_points: %d\n", min_points);
fprintf("  r_hss: %d\n", r_hss);
fprintf("  tol_hss: %.1e\n", tol_hss);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result = struct();

result.n = n;
result.N = N;
result.r_bf = r_bf;
result.tol_bf = tol_bf;
result.min_points = min_points;
result.r_hss = r_hss;
result.tol_hss = tol_hss;
result.num_sample = num_sample;
result.tol_cg = tol_cg;
result.maxit_cg = maxit_cg;

half_n = n / 2;

% Initialization.
x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);
xi_co = (-half_n : (half_n - 1))';
xi = TensorProduct2D(xi_co, xi_co);

% BF.
fprintf("BF.\n");
tic;
[K_BF, ~] = fastBF(exp_phi_func, x, xi, r_bf, tol_bf);
result.t_BF_construct = toc;
fprintf("  t_BF_construct: %.1e\n", result.t_BF_construct);

f_ex = randn(N, 1) + 1i * randn(N, 1);
tic;
Kf = apply_fbf(K_BF, f_ex);
result.t_BF_apply = toc;

fprintf("  t_BF_apply: %.1e\n", result.t_BF_apply);
result.rel_err_BF = fbf_check(N, exp_phi_func, f_ex, x, xi, Kf, num_sample);
fprintf("  rel_err_BF: %.1e\n", result.rel_err_BF);

% HSS.
fprintf("HSS.\n");
op_G = @(v) apply_fbf_adj(K_BF, apply_fbf(K_BF, v));
Gf_ex = op_G(f_ex);

result.rel_err_HSS = inf;
while result.rel_err_HSS >= tol_hss * 10
    fprintf("  current r_hss: %d\n", r_hss);
    tic;
    G_HSS = BF_HSS2D(n, n);
    G_HSS.BuildTree(min_points);
    G_HSS.BlackBoxConstruct_FastBF(K_BF, r_hss, tol_hss);
    result.t_HSS_construct = toc;
    fprintf("  t_HSS_construct: %.1e\n", result.t_HSS_construct);
    
    Gf = G_HSS.MyApply(f_ex);
    result.rel_err_HSS = norm(Gf - Gf_ex) / norm(Gf);
    fprintf("  rel_err_HSS: %.1e\n", result.rel_err_HSS);
    
    result.hss_rank = G_HSS.Rank();
    fprintf("  HSS rank: %d\n", result.hss_rank);
    result.hss_mem = G_HSS.Storage();
    ratio = result.hss_mem / N^2;
    fprintf("  ratio: %.1e\n", ratio);
    r_hss = r_hss * 2;
end

tic;
G_HSS.ULV_Factor();
result.t_HSS_factor = toc;
fprintf("  t_HSS_factor: %.1e\n", result.t_HSS_factor);

rhs = apply_fbf_adj(K_BF, Kf);
% Direct Solution.
fprintf("Direct solution.\n");
tic;
f_direct = G_HSS.Solve(rhs);
result.t_solve_direct = toc;

fprintf("  t_solve_direct: %.1e\n", result.t_solve_direct);
result.rel_res_direct = norm(Kf - apply_fbf(K_BF, f_direct)) / norm(Kf);
fprintf("  rel_res_direct: %.1e\n", result.rel_res_direct);
result.rel_err_direct = norm(f_ex - f_direct) / norm(f_ex);
fprintf("  rel_err_direct: %.1e\n", result.rel_err_direct);

% Iterative solution.
fprintf("Iterative solution.\n");
fprintf("  Without precond.\n");
tic;
[f_cg, result.flag_cg, ~, result.iter_cg] = pcg(op_G, rhs, tol_cg, maxit_cg);
result.t_cg = toc;

fprintf("    t_cg: %.1e\n", result.t_cg);
fprintf("    iter_cg: %d\n", result.iter_cg);
result.rel_res_cg = norm(Kf - apply_fbf(K_BF, f_cg)) / norm(Kf);
fprintf("    rel_res_cg: %.1e\n", result.rel_res_cg);
result.rel_err_cg = norm(f_ex - f_cg) / norm(f_ex);
fprintf("    rel_err_cg: %.1e\n", result.rel_err_cg);

fprintf("  With precond.\n");
op_M = @(v) G_HSS.Solve(v);
tic;
[f_pcg, result.flag_pcg, ~, result.iter_pcg] = pcg(op_G, rhs, tol_cg, maxit_cg, op_M);
result.t_pcg = toc;

fprintf("    t_pcg: %.1e\n", result.t_pcg);
fprintf("    iter_pcg: %d\n", result.iter_pcg);
result.rel_res_pcg = norm(Kf - apply_fbf(K_BF, f_pcg)) / norm(Kf);
fprintf("    rel_res_pcg: %.1e\n", result.rel_res_pcg);
result.rel_err_pcg = norm(f_ex - f_pcg) / norm(f_ex);
fprintf("    rel_err_pcg: %.1e\n", result.rel_err_pcg);

end