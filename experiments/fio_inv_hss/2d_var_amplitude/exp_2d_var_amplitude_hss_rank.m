exp_2d_var_amplitude_settings;

for it_p = 1 : num_n
    p = p_list(it_p);

    n = 2^p;
    N = n^2;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = run_FIO2D_hss_rank(...
        k_func, n, ...
        rank_func_tol);
    save("./data/2d_var_amplitude_hss_rank" ...
        + "_" + string(p) ...
        + ".mat", "result");
end