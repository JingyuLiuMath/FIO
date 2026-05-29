function result = run_FIO1D_var_hss_rank(...
    m, sigma_sq, ...
    exp_phi_func, N)

fprintf("Basic info.\n");
fprintf("  N: %d\n", N);

fprintf("  m: %d\n", m);
fprintf("  sigma_sq: %.1e\n", sigma_sq);

result = struct();

result.N = N;
result.tol_list = [1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8];
result.rank_list = zeros(size(result.tol_list));
result.rank_est_list = zeros(size(result.tol_list));
half_N = N / 2;

x_pts = rand(m, 1);
xi_pts = (rand(m, 1) - 0.5);

b_func = @(x) exp(-(x - x_pts.').^2 / sigma_sq);
c_func = @(xi) exp(-(xi / N - xi_pts.').^2 / sigma_sq);
a_func = @(x, xi) b_func(x) * c_func(xi).';

% Initialization.
x = (0 : (N - 1))' / N;

xi_row = (-half_N : -1)';
xi_col = (0 : (half_N - 1))';
sub_K = (exp_phi_func(x, xi_row) .* a_func(x, xi_row))' ...
    * (exp_phi_func(x, xi_col) .* a_func(x, xi_col));
sigma = svd(sub_K);
for it = 1 : size(result.tol_list, 1)
    tol = result.tol_list(it);
    result.rank_list(it) = find(sigma >= tol * sigma(1), 1, "last");
    result.rank_est_list(it) = ceil(log10(1 / tol) * log2(N) / 2);
    fprintf("  tol: %.1e, " + ...
        "rank: %d, " + ...
        "rank_est: %d\n", ...
        tol, ...
        result.rank_list(it), ...
        result.rank_est_list(it));
end


end