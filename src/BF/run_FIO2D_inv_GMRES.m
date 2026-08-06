function result = run_FIO2D_inv_GMRES(...
    n, ...
    a_func, phi_func, ...
    n_leaf_bf, r_bf, tol_bf, type_bf, ...
    num_sample, ...
    restart_gmres, tol_gmres, maxit_gmres)

N = n^2;
type_bf = string(type_bf);

fprintf("Basic info.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);

fprintf("  n_leaf_bf: %d\n", n_leaf_bf);
fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);
fprintf("  type_bf: %s\n", type_bf);

fprintf("  restart_gmres: %d\n", restart_gmres);
fprintf("  tol_gmres: %.1e\n", tol_gmres);
fprintf("  maxit_gmres: %d\n", maxit_gmres);

result = struct();

result.n = n;
result.N = N;
result.n_leaf_bf = n_leaf_bf;
result.r_bf = r_bf;
result.tol_bf = tol_bf;
result.type_bf = type_bf;
result.num_sample = num_sample;
result.restart_gmres = restart_gmres;
result.tol_gmres = tol_gmres;
result.maxit_gmres = maxit_gmres;

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

% FastBF.
fprintf("BF.\n");
t_construct_BF_start = tic;
result.K_BF = my_construct_bf(...
    n, a_func, phi_func, ...
    x, xi, ...
    n_leaf_bf, r_bf, tol_bf, type_bf);
result.t_construct_BF = toc(t_construct_BF_start);
fprintf("  t_construct_BF: %.1e\n", result.t_construct_BF);

result.ratio_BF = my_nnz_bf(result.K_BF, type_bf) / N^2;
fprintf("  ratio_BF: %.1e\n", result.ratio_BF);

result.f_ex = randn(N, 1) + 1i * randn(N, 1);
t_apply_BF_start = tic;
result.Kf_ex = my_apply_bf(...
    result.K_BF, result.f_ex, result.type_bf);
result.t_apply_BF = toc(t_apply_BF_start);

fprintf("  t_apply_BF: %.1e\n", result.t_apply_BF);
result.rel_err_BF = my_check_bf(...
    N, k_func, result.f_ex, x, xi, result.Kf_ex, num_sample);
fprintf("  rel_err_BF: %.1e\n", result.rel_err_BF);

% Iterative solution.
fprintf("Iterative solution.\n");
fprintf("  GMRES without precond.\n");
op_K = @(v) my_apply_bf(...
    result.K_BF, v, result.type_bf);
t_gmres_start = tic;
[f_gmres, result.flag_gmres, ~, result.iter_gmres] = gmres(...
    op_K, result.Kf_ex, restart_gmres, tol_gmres, maxit_gmres);
result.t_gmres = toc(t_gmres_start);
result.num_iter_gmres = ...
    max(result.iter_gmres(1) - 1, 0) * restart_gmres ...
    + result.iter_gmres(2);

fprintf("    t_gmres: %.1e\n", result.t_gmres);
fprintf("    iter_gmres: %d\n", result.num_iter_gmres);
result.rel_res_gmres = norm(...
    result.Kf_ex - op_K(f_gmres)) ...
    / norm(result.Kf_ex);
fprintf("    rel_res_gmres: %.1e\n", result.rel_res_gmres);
result.rel_err_gmres = norm(result.f_ex - f_gmres) ...
    / norm(result.f_ex);
fprintf("    rel_err_gmres: %.1e\n", result.rel_err_gmres);

end
