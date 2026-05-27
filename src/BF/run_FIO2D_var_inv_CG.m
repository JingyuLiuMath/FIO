function result = run_FIO2D_var_inv_CG(...
    a_func, r_a, tol_a, ...
    exp_phi_func, ...
    n, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg)

N = n^2;
fprintf("Basic info.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);

fprintf("  r_a: %d\n", r_a);
fprintf("  tol_a: %.1e\n", tol_a);

fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);

fprintf("  tol_cg: %.1e\n", tol_cg);
fprintf("  maxit_cg: %d\n", maxit_cg);

result = struct();

result.n = n;
result.N = N;
result.r_a = r_a;
result.tol_a = tol_a;
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

result.f_ex = randn(N, 1) + 1i * randn(N, 1);

% Amplitude term.
fprintf("Low-rank factorization on amplitude.\n");
[result.B, result.C] = Construct_Amplitude_LR_2D(a_func, x, xi, r_a, tol_a);

result.Kf_ex = result.B * (result.C' * result.f_ex);

result.rank_amplitude = size(result.B, 2);
result.rel_err_amplitude = fbf_check(N, a_func, result.f_ex, x, xi, result.Kf_ex, num_sample);
fprintf("  rank_amplitude: %d\n", result.rank_amplitude);
fprintf("  rel_err_amplitude: %.1e\n", result.rel_err_amplitude);

% BF.
fprintf("BF.\n");
tic;
[result.K_BF, ~] = fastBF(exp_phi_func, x, xi, r_bf, tol_bf);
result.t_construct_BF = toc;
fprintf("  t_construct_BF: %.1e\n", result.t_construct_BF);

tic;
result.Kf_ex = apply_fbf(result.K_BF, result.f_ex);
result.t_apply_BF = toc;

fprintf("  t_apply_BF: %.1e\n", result.t_apply_BF);
result.rel_err_BF = fbf_check(N, exp_phi_func, result.f_ex, x, xi, result.Kf_ex, num_sample);
fprintf("  rel_err_BF: %.1e\n", result.rel_err_BF);

% Apply.
fprintf("Apply.\n");

apply_func = @(f) apply(result.B, result.K_BF, result.C, f);
apply_adj_func = @(f) apply_adj(result.B, result.K_BF, result.C, f);

result.Kf_ex = apply_func(result.f_ex);

k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);
result.rel_err_apply = fbf_check(N, k_func, result.f_ex, x, xi, result.Kf_ex, num_sample);
fprintf("  rel_err_apply: %.1e\n", result.rel_err_apply);

op_G = @(v) apply_adj_func(apply_func(v));
rhs = apply_adj_func(result.Kf_ex);

% Iterative solution.
fprintf("Iterative solution.\n");
fprintf("  Without precond.\n");
tic;
[f_cg, result.flag_cg, ~, result.iter_cg] = pcg(op_G, rhs, tol_cg, maxit_cg);
result.t_cg = toc;

fprintf("    t_cg: %.1e\n", result.t_cg);
fprintf("    iter_cg: %d\n", result.iter_cg);
result.rel_res_cg = norm(result.Kf_ex - apply_func(f_cg)) / norm(result.Kf_ex);
fprintf("    rel_res_cg: %.1e\n", result.rel_res_cg);
result.rel_err_cg = norm(result.f_ex - f_cg) / norm(result.f_ex);
fprintf("    rel_err_cg: %.1e\n", result.rel_err_cg);

end

function u = apply(B, K_BF, C, f)

m = size(C, 2);
u = 0;
for j = 1 : m
    u = u + B(:, j) .* apply_fbf_batch(K_BF, conj(C(:, j)) .* f);
end

end

function u = apply_adj(B, K_BF, C, f)

m = size(C, 2);
u = 0;
for j = 1 : m
    u = u + C(:, j) .* apply_fbf_adj_batch(K_BF, conj(B(:, j)) .* f);
end

end