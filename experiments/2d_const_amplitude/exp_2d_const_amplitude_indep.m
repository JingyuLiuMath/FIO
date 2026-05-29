exp_2d_const_amplitude_settings;

for it_tol_hss = 1 : num_tol_hss
    for it_p = num_n : -1 : 1
        tol_hss = tol_hss_list(it_tol_hss);
        p = p_list(it_p);

        n = 2^p;
        N = n^2;

        fprintf("\n\n\n\n");
        fprintf("p: %d\n", p);

        load(data_path + "2d_const_amplituide_results_cg" ...
            + "_" + string(p) ...
            + ".mat");

        result = run_FIO2D_inv(...
            result, ...
            min_points, rank_func_tol, tol_hss);

        save(data_path + "2d_const_amplituide_results_indep" ...
            + "_" + string(p) ...
            + "_" + string(tol_hss) ...
            + ".mat", "result");
    end
end

path(originalPath);
