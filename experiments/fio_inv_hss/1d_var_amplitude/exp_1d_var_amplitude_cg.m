exp_1d_var_amplitude_settings;

for it_p = 1 : num_n
    p = p_list(it_p);

    N = 2^p;
    
    b_func = @(x) exp(-(x - x_pts.').^2 / sigma_sq);
    c_func = @(xi) exp(-(xi / N - xi_pts.').^2 / sigma_sq);
    a_func = @(x, xi) b_func(x) * c_func(xi).';

    k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);
    
    result = [];
    result = run_FIO1D_inv_HSS_CG(...
        N, ...
        a_func, phi_func, ...
        n_leaf_bf, r_bf, tol_bf, type_bf, ...
        num_sample, ...
        tol_cg, maxit_cg);

    save(data_path + "1d_var_amplituide_results_cg" ...
        + "_" + string(p) ...
        + ".mat", "result", '-v7.3');
end

path(originalPath);
