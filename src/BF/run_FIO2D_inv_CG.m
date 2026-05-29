function result = run_FIO2D_inv_CG(...
    k_func, n, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg)

N = n^2;
fprintf("Basic info.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);

fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result = struct();

result.n = n;
result.N = N;
result.r_bf = r_bf;
result.tol_bf = tol_bf;
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
[result.K_BF, ~] = fastBF(k_func, x, xi, r_bf, tol_bf);
result.t_construct_BF = toc;
fprintf("  t_construct_BF: %.1e\n", result.t_construct_BF);

result.f_ex = randn(N, 1) + 1i * randn(N, 1);
tic;
result.Kf_ex = apply_fbf(result.K_BF, result.f_ex);
result.t_apply_BF = toc;

fprintf("  t_apply_BF: %.1e\n", result.t_apply_BF);
result.rel_err_BF = fbf_check(N, k_func, result.f_ex, x, xi, result.Kf_ex, num_sample);
fprintf("  rel_err_BF: %.1e\n", result.rel_err_BF);

op_G = @(v) apply_fbf_adj(result.K_BF, apply_fbf(result.K_BF, v));
rhs = apply_fbf_adj(result.K_BF, result.Kf_ex);

% Iterative solution.
fprintf("Iterative solution.\n");
fprintf("  Without precond.\n");
tic;
[f_cg, result.flag_cg, ~, result.iter_cg] = pcg(op_G, rhs, tol_cg, maxit_cg);
result.t_cg = toc;

fprintf("    t_cg: %.1e\n", result.t_cg);
fprintf("    iter_cg: %d\n", result.iter_cg);
result.rel_res_cg = norm(result.Kf_ex - apply_fbf(result.K_BF, f_cg)) / norm(result.Kf_ex);
fprintf("    rel_res_cg: %.1e\n", result.rel_res_cg);
result.rel_err_cg = norm(result.f_ex - f_cg) / norm(result.f_ex);
fprintf("    rel_err_cg: %.1e\n", result.rel_err_cg);

end