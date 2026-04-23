function result = run_FIO_hss_rank(...
    exp_phi_func, N)

result = struct();

result.N = N;
result.tol_list = [1e-2; 1e-4; 1e-6; 1e-8];
result.rank_list = zeros(size(result.tol_list));
half_N = N / 2;

% Initialization.
x = (0 : (N - 1))' / N;

xi_row = (-half_N : -1)';
xi_col = (0 : (half_N - 1))';
sub_K = exp_phi_func(x, xi_row)' * exp_phi_func(x, xi_col);
sigma = svd(sub_K);
for it = 1 : size(result.tol_list, 1)
    result.rank_list(it) = find(sigma >= result.tol_list(it) * sigma(1), 1, "last");
    fprintf("tol: %.1e, rank: %d\n", result.tol_list(it), result.rank_list(it));
end


end