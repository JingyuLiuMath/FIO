function result = run_FIO2D_inv_ID_Full(...
    exp_phi_func, ...
    result, ...
    min_points, tol_hss)

n = result.n;
N = result.N;

r_bf = result.r_bf;
tol_bf = result.tol_bf;

tol_cg = result.tol_cg;
maxit_cg = result.maxit_cg;

fprintf("Basic info.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);

fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);

fprintf("  min_points: %d\n", min_points);
fprintf("  tol_hss: %.1e\n", tol_hss);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result.min_points = min_points;
result.tol_hss = tol_hss;

K_BF = result.K_BF;
f_ex = result.f_ex;
Kf_ex = result.Kf_ex;
result = rmfield(result, {'K_BF', 'f_ex', 'Kf_ex'});

half_n = n / 2;

% Initialization.
x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);
xi_co = (-half_n : (half_n - 1))';
xi = TensorProduct2D(xi_co, xi_co);

op_G = @(v) apply_fbf_adj(K_BF, apply_fbf(K_BF, v));
rhs = apply_fbf_adj(K_BF, Kf_ex);

% HSS.
fprintf("HSS.\n");
Gf_ex = op_G(f_ex);

tic;
G_HSS = BF_HSS2D(n, n);
G_HSS.BuildTree(min_points);
G_HSS.Construct_ID_Full(exp_phi_func, x, xi(G_HSS.perm_, :), tol_hss);
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

tic;
G_HSS.ULV_Factor();
result.t_HSS_factor = toc;
fprintf("  t_HSS_factor: %.1e\n", result.t_HSS_factor);

% Direct Solution.
fprintf("Direct solution.\n");
tic;
f_direct = G_HSS.Solve(rhs);
result.t_solve_direct = toc;

fprintf("  t_solve_direct: %.1e\n", result.t_solve_direct);
result.rel_res_direct = norm(Kf_ex - apply_fbf(K_BF, f_direct)) / norm(Kf_ex);
fprintf("  rel_res_direct: %.1e\n", result.rel_res_direct);
result.rel_err_direct = norm(f_ex - f_direct) / norm(f_ex);
fprintf("  rel_err_direct: %.1e\n", result.rel_err_direct);

% Iterative solution.
fprintf("Iterative solution.\n");
fprintf("  With precond.\n");
op_M = @(v) G_HSS.Solve(v);
tic;
[f_pcg, result.flag_pcg, ~, result.iter_pcg] = pcg(op_G, rhs, tol_cg, maxit_cg, op_M);
result.t_pcg = toc;

fprintf("    t_pcg: %.1e\n", result.t_pcg);
fprintf("    iter_pcg: %d\n", result.iter_pcg);
result.rel_res_pcg = norm(Kf_ex - apply_fbf(K_BF, f_pcg)) / norm(Kf_ex);
fprintf("    rel_res_pcg: %.1e\n", result.rel_res_pcg);
result.rel_err_pcg = norm(f_ex - f_pcg) / norm(f_ex);
fprintf("    rel_err_pcg: %.1e\n", result.rel_err_pcg);

end