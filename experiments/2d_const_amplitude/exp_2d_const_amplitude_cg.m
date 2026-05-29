exp_2d_const_amplitude_settings;

for it_p = 1 : num_n
    p = p_list(it_p);

    n = 2^p;
    N = n^2;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = [];
    result = run_FIO2D_inv_CG(...
        k_func, n, ...
        r_bf, tol_bf, ...
        num_sample, ...
        tol_cg, maxit_cg);

    save(data_path + "2d_const_amplituide_results_cg" ...
        + "_" + string(p) ...
        + ".mat", "result", '-v7.3');
end

path(originalPath);
