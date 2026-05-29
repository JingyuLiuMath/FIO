exp_1d_var_amplitude_settings;

for it_p = 1 : num_n
    p = p_list(it_p);

    N = 2^p;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = [];
    result = run_FIO1D_var_inv_CG(...
        m, sigma_sq, ...
        exp_phi_func, N, ...
        r_bf, tol_bf, ...
        num_sample, ...
        tol_cg, maxit_cg);

    save(data_path + "1d_var_amplituide_results_cg" ...
        + "_" + string(p) ...
        + ".mat", "result", '-v7.3');
end

path(originalPath);
