function result = run_FIO1D_inv_NewtonSchulz(...
    result_bf, ...
    r_ns, tol_ns, maxit_ns, num_test, type_bf_inv)

N = result_bf.N;

n_leaf_bf = result_bf.n_leaf_bf;
r_bf = result_bf.r_bf;
tol_bf = result_bf.tol_bf;
type_bf = string(result_bf.type_bf);
type_bf_inv = string(type_bf_inv);

restart_gmres = result_bf.restart_gmres;
tol_gmres = result_bf.tol_gmres;
maxit_gmres = result_bf.maxit_gmres;

fprintf("Basic info.\n");
fprintf("  N: %d\n", N);

fprintf("  n_leaf_bf: %d\n", n_leaf_bf);
fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);
fprintf("  type_bf: %s\n", type_bf);

fprintf("  r_ns: %d\n", r_ns);
fprintf("  tol_ns: %.1e\n", tol_ns);
fprintf("  maxit_ns: %d\n", maxit_ns);
fprintf("  num_test: %d\n", num_test);
fprintf("  type_bf_inv: %s\n", type_bf_inv);

fprintf("  restart_gmres: %d\n", restart_gmres);
fprintf("  tol_gmres: %.1e\n", tol_gmres);
fprintf("  maxit_gmres: %d\n", maxit_gmres);

result = struct();

result.N = N;
result.n_leaf_bf = n_leaf_bf;
result.r_bf = r_bf;
result.tol_bf = tol_bf;
result.type_bf = type_bf;
result.r_ns = r_ns;
result.tol_ns = tol_ns;
result.maxit_ns = maxit_ns;
result.num_test = num_test;
result.type_bf_inv_requested = type_bf_inv;
result.type_bf_inv = type_bf_inv;
result.restart_gmres = restart_gmres;
result.tol_gmres = tol_gmres;
result.maxit_gmres = maxit_gmres;

result.t_construct_BF = result_bf.t_construct_BF;
result.t_apply_BF = result_bf.t_apply_BF;
result.ratio_BF = result_bf.ratio_BF;
result.rel_err_BF = result_bf.rel_err_BF;

result.t_gmres = result_bf.t_gmres;
result.flag_gmres = result_bf.flag_gmres;
result.iter_gmres = result_bf.iter_gmres;
result.num_iter_gmres = result_bf.num_iter_gmres;
result.rel_res_gmres = result_bf.rel_res_gmres;
result.rel_err_gmres = result_bf.rel_err_gmres;

op_K = @(v) my_apply_bf(...
    result_bf.K_BF, v, result.type_bf);

fprintf("  rel_err_BF: %.1e\n", result_bf.rel_err_BF);

% Newton-Schulz inverse.
fprintf("Newton-Schulz inverse.\n");
t_construct_NS_start = tic;
[result.K_inv_BF, result.history_ns] = ...
    my_bf_inv_newton_schulz(...
    N, result_bf.K_BF, n_leaf_bf, ...
    r_ns, tol_ns, maxit_ns, num_test, ...
    1, type_bf, type_bf_inv);
result.type_bf_inv = string(result.history_ns.output_type_bf);
result.t_construct_NS = toc(t_construct_NS_start);
fprintf("  t_construct_NS: %.1e\n", result.t_construct_NS);
fprintf("  num_iter_ns: %d\n", ...
    result.history_ns.num_iterations);
fprintf("  type_bf_inv_output: %s\n", result.type_bf_inv);

result.ratio_NS = result.history_ns.nnz(end) / N^2;
fprintf("  ratio_NS: %.1e\n", result.ratio_NS);

% Direct inverse application.
fprintf("Direct inverse application.\n");
t_apply_NS_start = tic;
f_ns = my_apply_bf(...
    result.K_inv_BF, result_bf.Kf_ex, result.type_bf_inv);
result.t_apply_NS = toc(t_apply_NS_start);
fprintf("  t_apply_NS: %.1e\n", result.t_apply_NS);

result.rel_res_ns = norm(...
    result_bf.Kf_ex - op_K(f_ns)) / norm(result_bf.Kf_ex);
fprintf("  rel_res_ns: %.1e\n", result.rel_res_ns);
result.rel_err_ns = norm(result_bf.f_ex - f_ns) ...
    / norm(result_bf.f_ex);
fprintf("  rel_err_ns: %.1e\n", result.rel_err_ns);

% Preconditioned iterative solution.
fprintf("Iterative solution.\n");
fprintf("  GMRES with Newton-Schulz precond.\n");
op_M = @(v) my_apply_bf(...
    result.K_inv_BF, v, result.type_bf_inv);
t_pgmres_start = tic;
[f_pgmres, result.flag_pgmres, ~, result.iter_pgmres] = gmres(...
    op_K, result_bf.Kf_ex, restart_gmres, ...
    tol_gmres, maxit_gmres, op_M);
result.t_pgmres = toc(t_pgmres_start);
result.num_iter_pgmres = ...
    max(result.iter_pgmres(1) - 1, 0) * restart_gmres ...
    + result.iter_pgmres(2);

fprintf("    t_pgmres: %.1e\n", result.t_pgmres);
fprintf("    iter_pgmres: %d\n", result.num_iter_pgmres);
result.rel_res_pgmres = norm(...
    result_bf.Kf_ex - op_K(f_pgmres)) / norm(result_bf.Kf_ex);
fprintf("    rel_res_pgmres: %.1e\n", result.rel_res_pgmres);
result.rel_err_pgmres = norm(result_bf.f_ex - f_pgmres) ...
    / norm(result_bf.f_ex);
fprintf("    rel_err_pgmres: %.1e\n", result.rel_err_pgmres);

end
