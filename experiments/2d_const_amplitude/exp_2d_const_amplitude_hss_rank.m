clear;
close all;
warning off;

addpath('../../extern/FastBF.m/src');
addpath('../../extern/FastBF.m/test/kernels');

exp_phi_func = @(x, xi) fun0_2D(x, xi);

p_list = (4 : 7)';
num_n = length(p_list);

for it_p = 1 : num_n
    p = p_list(it_p);

    n = 2^p;
    N = n^2;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = run_FIO_hss_rank2D(exp_phi_func, n);
    save("./data/hss_rank_" + string(p) + ".mat", ...
        "result");
end