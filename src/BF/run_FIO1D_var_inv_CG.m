function result = run_FIO1D_var_inv_CG(...
    m, sigma_sq, ...
    exp_phi_func, N, ...
    r_bf, tol_bf, ...
    num_sample, ...
    tol_cg, maxit_cg)

fprintf("Basic info.\n");
fprintf("  N: %d\n", N);

fprintf("  m: %d\n", m);
fprintf("  sigma_sq: %.1e\n", sigma_sq);

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

x_pts = rand(m, 1);
xi_pts = rand(m, 1) * (1 - 1 / N) - 0.5;
% x_pts = linspace(0, 1, m)';
% x_pts = x_pts + randn(size(x_pts)) * 1e-6;
% xi_pts = linspace(-0.5, 0.5 - 1 / N, m)';
% xi_pts = xi_pts + randn(size(xi_pts)) * 1e-6;

result.m = m;
result.sigma_sq = sigma_sq;
result.x_pts = x_pts;
result.xi_pts = xi_pts;

b_func = @(x) exp(-(x - x_pts.').^2 / sigma_sq);
c_func = @(xi) exp(-(xi / N - xi_pts.').^2 / sigma_sq);
a_func = @(x, xi) b_func(x) * c_func(xi).';

% Initialization.
x = (0 : (N - 1))' / N;
xi = (-half_N : (half_N - 1))';

result.f_ex = randn(N, 1) + 1i * randn(N, 1);

% Amplitude term.
fprintf("Low-rank factorization on amplitude.\n");
result.B = b_func(x);
result.C = conj(c_func(xi));

result.Kf_ex = result.B * (result.C' * result.f_ex);

result.rank_amplitude = m;
fprintf("  rank_amplitude: %d\n", result.rank_amplitude);

% BF.
fprintf("BF.\n");
t_construct_BF_start = tic;
[result.K_BF, ~] = fastBF(exp_phi_func, x, xi, r_bf, tol_bf);
result.t_construct_BF = toc(t_construct_BF_start);
fprintf("  t_construct_BF: %.1e\n", result.t_construct_BF);

M = size(result.K_BF.M, 1);
fprintf("  M in BF: %d\n", M);

t_apply_BF_start = tic;
result.Kf_ex = apply_fbf(result.K_BF, result.f_ex);
result.t_apply_BF = toc(t_apply_BF_start);

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
t_cg_start = tic;
[f_cg, result.flag_cg, ~, result.iter_cg] = pcg(op_G, rhs, tol_cg, maxit_cg);
result.t_cg = toc(t_cg_start);

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