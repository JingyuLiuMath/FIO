clear;
close all;
warning off;

addpath('../../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_2D(x, xi);

p_list = (6 : 9)';
num_n = length(p_list);

r_bf = 30;
tol_bf = 1e-10;

min_points = 256;
tol_hss_list = [1e-3];
num_tol_hss = length(tol_hss_list);

tol_cg = 1e-12;
maxit_cg = 500;

num_sample = 256;

for it_tol_hss = 1 : num_tol_hss
    for it_p = 1 : num_n
        tol_hss = tol_hss_list(it_tol_hss);
        p = p_list(it_p);
    
        n = 2^p;
        N = n^2;     
        r_hss = 10 * n;
    
        fprintf("\n\n\n\n");
        fprintf("p: %d\n", p);
    
        result = run_FIO_inv_fastBF2D(...
            exp_phi_func, n, ...
            r_bf, tol_bf, ...
            min_points, r_hss, tol_hss, ...
            num_sample, ...
            tol_cg, maxit_cg);
        
        save("./data/2d_const_amplituide_results_" + string(p) + "_" + string(tol_hss) + ".mat", ...
            "result");
    end
end