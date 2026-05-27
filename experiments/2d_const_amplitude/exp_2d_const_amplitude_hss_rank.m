exp_2d_const_amplitude_settings;

for it_p = 1 : num_n
    p = p_list(it_p);

    n = 2^p;
    N = n^2;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = run_FIO2D_hss_rank(exp_phi_func, n);
    save("./data/2d_const_amplitude_hss_rank" ...
        + "_" + string(p) ...
        + ".mat", "result");
end