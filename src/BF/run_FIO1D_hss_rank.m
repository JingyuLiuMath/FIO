function result = run_FIO1D_hss_rank(...
    k_func, N, ...
    rank_func_tol)

fprintf("Basic info.\n");
fprintf("  N: %d\n", N);

result = struct();

result.N = N;
result.tol_list = [1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8];
result.rank_list = zeros(size(result.tol_list));
result.rank_est_list = zeros(size(result.tol_list));
half_N = N / 2;

% Initialization.
x = (0 : (N - 1))' / N;

xi_row = (-half_N : -1)';
xi_col = (0 : (half_N - 1))';
sub_K = k_func(x, xi_row)' * k_func(x, xi_col);
sigma = svd(sub_K);
for it = 1 : size(result.tol_list, 1)
    tol = result.tol_list(it);
    result.rank_list(it) = find(sigma >= tol * sigma(1), 1, "last");
    result.rank_est_list(it) = rank_func_tol(tol);
    fprintf("  tol: %.1e, " + ...
        "rank: %d, " + ...
        "rank_est: %d\n", ...
        tol, ...
        result.rank_list(it), ...
        result.rank_est_list(it));
end


end