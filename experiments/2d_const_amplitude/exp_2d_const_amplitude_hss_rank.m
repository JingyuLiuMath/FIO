clear;
close all;

exp_phi_func = @(x, xi) fun_2D(x, xi);

p_list = (4 : 7)';
num_n = length(p_list);

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