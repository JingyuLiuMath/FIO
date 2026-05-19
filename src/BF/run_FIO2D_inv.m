function result = run_FIO2D_inv(...
    result_bf, ...
    min_points, tol_hss, indep)

n = result_bf.n;
N = result_bf.N;

r_bf = result_bf.r_bf;
tol_bf = result_bf.tol_bf;

tol_cg = result_bf.tol_cg;
maxit_cg = result_bf.maxit_cg;

fprintf("Basic info.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);

fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);

fprintf("  min_points: %d\n", min_points);
fprintf("  tol_hss: %.1e\n", tol_hss);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

fprintf("  indep: %d\n", indep);

result.n = n;
result.N = N;
result.r_bf = r_bf;
result.tol_bf = tol_bf;
result.min_points = min_points;
result.tol_hss = tol_hss;
result.tol_cg = tol_cg;
result.maxit_cg = maxit_cg;
result.indep = indep;

op_G = @(v) apply_fbf_adj_batch(result_bf.K_BF, apply_fbf_batch(result_bf.K_BF, v));
rhs = apply_fbf_adj(result_bf.K_BF, result_bf.Kf_ex);

M = size(result_bf.K_BF.M, 1);
fprintf("  M in BF: %d\n", M);
fprintf("  rel_err_BF: %.1e\n", result_bf.rel_err_BF);

% HSS.
fprintf("HSS.\n");
Gf_ex = op_G(result_bf.f_ex);

G_HSS = BF_HSS2D(n, n);
G_HSS.BuildTree(min_points);

c = ceil(3 * 2 * log10(1 / tol_hss));
rank_func = @(ell) c * n / 2^ell;
result.rel_err_HSS = inf;
while result.rel_err_HSS >= tol_hss * 10
    fprintf("  c: %d\n", c);
    t_HSS_construct_start = tic;
    if indep == 1
        G_HSS.BlackBoxConstruct_Indep_FastBF(op_G, rank_func, tol_hss);
    else
        G_HSS.BlackBoxConstruct_FastBF(op_G, rank_func, tol_hss);
    end
    result.t_HSS_construct = toc(t_HSS_construct_start);
    fprintf("  t_HSS_construct: %.1e\n", result.t_HSS_construct);

    Gf = G_HSS.MyApply(result_bf.f_ex);
    result.rel_err_HSS = norm(Gf - Gf_ex) / norm(Gf);
    fprintf("  rel_err_HSS: %.1e\n", result.rel_err_HSS);

    result.hss_rank = G_HSS.Rank();
    fprintf("  HSS rank: %d\n", result.hss_rank);
    result.hss_mem = G_HSS.Storage();
    ratio = result.hss_mem / N^2;
    fprintf("  ratio: %.1e\n", ratio);

    c = c * 2;
    rank_func = @(ell) c * n / 2^ell;
end

t_HSS_factor_start = tic;
G_HSS.ULV_Factor();
result.t_HSS_factor = toc(t_HSS_factor_start);
fprintf("  t_HSS_factor: %.1e\n", result.t_HSS_factor);

% Direct Solution.
fprintf("Direct solution.\n");
t_direct_start = tic;
f_direct = G_HSS.Solve(rhs);
result.t_direct = toc(t_direct_start);

fprintf("  t_direct: %.1e\n", result.t_direct);
result.rel_res_direct = norm(result_bf.Kf_ex - apply_fbf(result_bf.K_BF, f_direct)) / norm(result_bf.Kf_ex);
fprintf("  rel_res_direct: %.1e\n", result.rel_res_direct);
result.rel_err_direct = norm(result_bf.f_ex - f_direct) / norm(result_bf.f_ex);
fprintf("  rel_err_direct: %.1e\n", result.rel_err_direct);

% Iterative solution.
fprintf("Iterative solution.\n");
fprintf("  With precond.\n");
op_M = @(v) G_HSS.Solve(v);
t_pcg_start = tic;
[f_pcg, result.flag_pcg, ~, result.iter_pcg] = pcg(op_G, rhs, tol_cg, maxit_cg, op_M);
result.t_pcg = toc(t_pcg_start);

fprintf("    t_pcg: %.1e\n", result.t_pcg);
fprintf("    iter_pcg: %d\n", result.iter_pcg);
result.rel_res_pcg = norm(result_bf.Kf_ex - apply_fbf(result_bf.K_BF, f_pcg)) / norm(result_bf.Kf_ex);
fprintf("    rel_res_pcg: %.1e\n", result.rel_res_pcg);
result.rel_err_pcg = norm(result_bf.f_ex - f_pcg) / norm(result_bf.f_ex);
fprintf("    rel_err_pcg: %.1e\n", result.rel_err_pcg);

end