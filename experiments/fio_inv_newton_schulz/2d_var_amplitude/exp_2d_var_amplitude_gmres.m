exp_2d_var_amplitude_settings;

for it_p = 1 : num_n
    p = p_list(it_p);
    n = 2^p;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = run_FIO2D_inv_GMRES( ...
        n, a_func, phi_func, ...
        n_leaf_bf, r_bf, tol_bf, type_bf, num_sample, ...
        restart_gmres, tol_gmres, maxit_gmres);

    save(data_path + "2d_var_amplitude_results_gmres" ...
        + "_" + string(p) + ".mat", "result", "-v7.3");
end
