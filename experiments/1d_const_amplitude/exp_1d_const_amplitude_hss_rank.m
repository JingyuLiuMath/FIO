exp_1d_const_amplitude_settings;

p_list = (10 : 12)';
num_n = length(p_list);

for it_p = 1 : num_n
    p = p_list(it_p);

    N = 2^p;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);
    
    result = run_FIO1D_hss_rank(...
        k_func, N, ...
        rank_func_tol);
    save("./data/1d_const_amplitude_hss_rank" ...
        + "_" + string(p) ...
        + ".mat", "result");
end

path(originalPath);
