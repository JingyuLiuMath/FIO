exp_1d_var_amplitude_settings;

p_list = (8 : 12)';
num_n = length(p_list);

for it_p = 1 : num_n
    p = p_list(it_p);

    N = 2^p;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = run_FIO1D_var_hss_rank(...
        m, sigma_sq, ...
        exp_phi_func, N);
    save("./data/1d_var_amplitude_hss_rank" ...
        + "_" + string(p) ...
        + ".mat", "result");
end

path(originalPath);
