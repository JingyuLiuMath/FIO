function result = run_FIO1D_inv_CG(...
    exp_phi_func, N, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg)

fprintf("Basic info.\n");
fprintf("  N: %d\n", N);

fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result = struct();

result.N = N;
result.r_bf = r_bf;
result.tol_bf = tol_bf;
result.num_sample = num_sample;
result.tol_cg = tol_cg;
result.maxit_cg = maxit_cg;

half_N = N / 2;

% Initialization.
x = (0 : (N - 1))' / N;
xi = (-half_N : (half_N - 1))';

% BF.
fprintf("BF.\n");
t_construct_BF_start = tic;
[result.K_BF, ~] = fastBF(exp_phi_func, x, xi, r_bf, tol_bf);
result.t_construct_BF = toc(t_construct_BF_start);
fprintf("  t_construct_BF: %.1e\n", result.t_construct_BF);

M = size(result.K_BF.M, 1);
fprintf("  M in BF: %d\n", M);

result.f_ex = randn(N, 1) + 1i * randn(N, 1);
t_apply_BF_start = tic;
result.Kf_ex = apply_fbf(result.K_BF, result.f_ex);
result.t_apply_BF = toc(t_apply_BF_start);

fprintf("  t_apply_BF: %.1e\n", result.t_apply_BF);
result.rel_err_BF = fbf_check(N, exp_phi_func, result.f_ex, x, xi, result.Kf_ex, num_sample);
fprintf("  rel_err_BF: %.1e\n", result.rel_err_BF);

op_G = @(v) apply_fbf_adj(result.K_BF, apply_fbf(result.K_BF, v));
rhs = apply_fbf_adj(result.K_BF, result.Kf_ex);

% Iterative solution.
fprintf("Iterative solution.\n");
fprintf("  Without precond.\n");
t_cg_start = tic;
[f_cg, result.flag_cg, ~, result.iter_cg] = pcg(op_G, rhs, tol_cg, maxit_cg);
result.t_cg = toc(t_cg_start);

fprintf("    t_cg: %.1e\n", result.t_cg);
fprintf("    iter_cg: %d\n", result.iter_cg);
result.rel_res_cg = norm(result.Kf_ex - apply_fbf(result.K_BF, f_cg)) / norm(result.Kf_ex);
fprintf("    rel_res_cg: %.1e\n", result.rel_res_cg);
result.rel_err_cg = norm(result.f_ex - f_cg) / norm(result.f_ex);
fprintf("    rel_err_cg: %.1e\n", result.rel_err_cg);

end