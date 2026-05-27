clear;
close all;

exp_phi_func = @(x, xi) fun_2D_2(x, xi);
a_func = @(x, xi) a_fun_2D_2(x, xi);

p_list = (4 : 6)';
num_n = length(p_list);

for it_p = 1 : num_n
    p = p_list(it_p);

    n = 2^p;
    N = n^2;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = run_FIO2D_var_hss_rank(exp_phi_func, n, a_func);
    save("./data/2d_var_amplitude_hss_rank" ...
        + "_" + string(p) ...
        + ".mat", "result");
end