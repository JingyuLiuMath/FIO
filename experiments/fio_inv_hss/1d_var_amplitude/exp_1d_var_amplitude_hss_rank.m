exp_1d_var_amplitude_settings;

p_list = (10 : 12)';
num_n = length(p_list);

for it_p = 1 : num_n
    p = p_list(it_p);

    N = 2^p;

    b_func = @(x) exp(-(x - x_pts.').^2 / sigma_sq);
    c_func = @(xi) exp(-(xi / N - xi_pts.').^2 / sigma_sq);
    a_func = @(x, xi) b_func(x) * c_func(xi).';

    k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);
    
    result = run_FIO1D_hss_rank(...
        k_func, N, ...
        rank_func_tol);
    save("./data/1d_var_amplitude_hss_rank" ...
        + "_" + string(p) ...
        + ".mat", "result");
end

path(originalPath);
