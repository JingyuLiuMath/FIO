function result = run_FIO2D_inv_HSS_CG(...
    n, ...
    a_func, phi_func, ...
    n_leaf_bf, r_bf, tol_bf, type_bf, ...
    num_sample, ...
    tol_cg, maxit_cg)

N = n^2;
fprintf("Basic info.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);

fprintf("  n_leaf_bf: %d\n", n_leaf_bf);
fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);
fprintf("  type_bf: %s\n", type_bf);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result = struct();

result.n = n;
result.N = N;
result.n_leaf_bf = n_leaf_bf;
result.r_bf = r_bf;
result.tol_bf = tol_bf;
result.type_bf = type_bf;
result.num_sample = num_sample;
result.tol_cg = tol_cg;
result.maxit_cg = maxit_cg;

half_n = n / 2;

% Initialization.
x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);
xi_co = (-half_n : (half_n - 1))';
xi = TensorProduct2D(xi_co, xi_co);

exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));
k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

% BF.
fprintf("BF.\n");
tic;
result.K_BF = my_construct_bf(...
    n, a_func, phi_func, ...
    x, xi, ...
    n_leaf_bf, r_bf, tol_bf, type_bf);
result.t_construct_BF = toc;
fprintf("  t_construct_BF: %.1e\n", result.t_construct_BF);

result.f_ex = randn(N, 1) + 1i * randn(N, 1);
tic;
result.Kf_ex = my_apply_bf(result.K_BF, result.f_ex, type_bf);
result.t_apply_BF = toc;

fprintf("  t_apply_BF: %.1e\n", result.t_apply_BF);
result.rel_err_BF = fbf_check(N, k_func, result.f_ex, x, xi, result.Kf_ex, num_sample);
fprintf("  rel_err_BF: %.1e\n", result.rel_err_BF);

op_G = @(v) my_apply_bf_adj(result.K_BF, ...
    my_apply_bf(result.K_BF, v, type_bf), type_bf);
rhs = my_apply_bf_adj(result.K_BF, result.Kf_ex, type_bf);

% Iterative solution.
fprintf("Iterative solution.\n");
fprintf("  Without precond.\n");
tic;
[f_cg, result.flag_cg, ~, result.iter_cg] = pcg(op_G, rhs, tol_cg, maxit_cg);
result.t_cg = toc;

fprintf("    t_cg: %.1e\n", result.t_cg);
fprintf("    iter_cg: %d\n", result.iter_cg);
result.rel_res_cg = norm(result.Kf_ex - my_apply_bf(result.K_BF, f_cg, type_bf)) / norm(result.Kf_ex);
fprintf("    rel_res_cg: %.1e\n", result.rel_res_cg);
result.rel_err_cg = norm(result.f_ex - f_cg) / norm(result.f_ex);
fprintf("    rel_err_cg: %.1e\n", result.rel_err_cg);

end
