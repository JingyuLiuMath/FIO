clear;
close all;
warning off;

originalPath = path;
addpath('../../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_1D(x, xi);

p_list = (7 : 14)';
num_n = length(p_list);

for it_p = 1 : num_n
    p = p_list(it_p);

    N = 2^p;

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result = run_FIO1D_hss_rank(exp_phi_func, N);
    save("./data/hss_rank" ...
        + "_" + string(p) ...
        + ".mat", "result");
end

path(originalPath);
