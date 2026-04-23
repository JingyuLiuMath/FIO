clear;
close all;
warning off;

addpath('../../extern/FastBF.m/src');
addpath('../../extern/FastBF.m/test/kernels');

exp_phi_func = @(x, xi) fun0_1D(x, xi);

p_list = (10 : 2 : 18)';
% p_list = (10 : 2 : 12)';
num_n = length(p_list);

r_bf = 10;
tol_bf = 1e-8;

min_points = 256;
tol_hss_list = [1e-2, 1e-4];
num_tol_hss = length(tol_hss_list);

tol_cg = 1e-12;
maxit_cg = 500;

num_sample = 256;

warmup = 1;
for it_tol_hss = 1 : num_tol_hss
    for it_p = 1 : num_n
        tol_hss = tol_hss_list(it_tol_hss);
        p = p_list(it_p);

        N = 2^p;
        r_hss = 2 * p;

        fprintf("\n\n\n\n");
        fprintf("Basic info.\n");
        fprintf("  p: %d\n", p);
        fprintf("  N: %d\n", N);

        fprintf("  r_bf: %d\n", r_bf);
        fprintf("  tol_bf: %.1e\n", tol_bf);

        fprintf("  min_points: %d\n", min_points);
        fprintf("  r_hss: %d\n", r_hss);
        fprintf("  tol_hss: %.1e\n", tol_hss);

        fprintf("  tol_cg: %.1e\n", tol_cg);
        fprintf("  maxit_cg: %d\n", maxit_cg);
        result = run_FIO_inv_fastBF(...
            exp_phi_func, N, ...
            r_bf, tol_bf, ...
            min_points, r_hss, tol_hss, ...
            num_sample, ...
            tol_cg, maxit_cg);

        save("./data/1d_const_amplituide_results_" + string(p) + "_" + string(tol_hss) + ".mat", ...
            "result");
    end
end