function result = run_FIO2D_hss_rank(...
    exp_phi_func, n)

N = n^2;
fprintf("Basic info.\n");
fprintf("  n: %d\n", n);

result = struct();

result.n = n;
result.N = N;
result.tol_list = [1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8];
result.rank_list = zeros(size(result.tol_list));
result.rank_est_list = zeros(size(result.tol_list));
half_n = n / 2;

% Initialization.
x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);

xi_co1 = (-half_n : -1)';
xi_co2 = (0 : (half_n - 1))';

xi_row = TensorProduct2D(xi_co1, xi_co1);
xi_col = TensorProduct2D(xi_co1, xi_co2);

factor = 3;
sub_K = exp_phi_func(x, xi_row)' * exp_phi_func(x, xi_col);
sigma = svd(sub_K);
for it = 1 : size(result.tol_list, 1)
    tol = result.tol_list(it);
    result.rank_list(it) = find(sigma >= tol * sigma(1), 1, "last");
    result.rank_est_list(it) = ceil(factor * log10(1 / tol)) * half_n;
    fprintf("tol: %.1e, " + ...
        "rank: %d, " + ...
        "rank_est: %d\n", ...
        result.tol_list(it), ...
        result.rank_list(it), ...
        result.rank_est_list(it));
end


end