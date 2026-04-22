clear;
close all;
warning off;

c_func = @(x) (2 + sin(2 * pi * x)) / 8;
phi_func = @(x, xi) x * xi.' + c_func(x) * abs(xi).';

% p_list = (10 : 2 : 18)';
p_list = (10 : 2 : 12)';
num_n = length(p_list);

r_bf = 8;
tol_bf = 1e-12;

min_points = 256;
tol_hss_list = [1e-3, 1e-6];
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
        fprintf("Basic info.\n")
        fprintf("  p: %d\n", p);
        fprintf("  N: %d\n", N);
        
        fprintf("  r_bf: %d\n", r_bf);
        fprintf("  tol_bf: %.1e\n", tol_bf);

        fprintf("  min_points: %d\n", min_points);
        fprintf("  r_hss: %d\n", r_hss);
        fprintf("  tol_hss: %.1e\n", tol_hss);

        fprintf("  tol_cg: %.1e\n", tol_cg);
        fprintf("  maxit_cg: %d\n", maxit_cg);

        % BF.
        BF_input = struct();
        BF_input.xbox = [1,N+1];
        BF_input.xx = (1:N)';
        BF_input.kbox = [1,N+1];
        BF_input.kk = (1:N)';
        BF_input.fun = @(x,k)fun0(N,x,k);
        BF_input.mR = r_bf;
        BF_input.tol = tol_bf;

        % BF.
        fprintf("BF.\n");
        if warmup == 1
            K_BF = bf_explicit(BF_input.fun, ...
                BF_input.xx, BF_input.xbox, ...
                BF_input.kk, BF_input.kbox, ...
                BF_input.mR, BF_input.tol, 0);
            warmup = 0;
        end

        tic;
        K_BF = bf_explicit(BF_input.fun, ...
            BF_input.xx, BF_input.xbox, ...
            BF_input.kk, BF_input.kbox, ...
            BF_input.mR, BF_input.tol, 0);
        t_BF_construct = toc;
        fprintf("  t_BF_construct: %.1e\n", t_BF_construct);

        f_ex = randn(N,1) + 1i*randn(N,1);
        tic;
        Kf = apply_bf(K_BF, f_ex);
        t_BF_apply = toc;

        fprintf("  t_BF_apply: %.1e\n", t_BF_apply);
        rel_err_BF = bf_explicit_check(N, BF_input.fun, f_ex, ...
            BF_input.xx, BF_input.kk, Kf, num_sample);
        fprintf("  rel_err_BF: %.1e\n", rel_err_BF);

        % HSS.
        fprintf("HSS.\n");
        op_G = @(v) apply_bf_adj(K_BF, apply_bf(K_BF, v));

        tic;
        G_HSS = BF_HSS(N);
        G_HSS.BuildTree(min_points);
        leaf_size = G_HSS.MaxLeafSize();
        G_HSS.BlackBoxConstruct(op_G, leaf_size, r_hss, tol_hss);
        t_HSS_construct = toc;
        fprintf("  t_HSS_construct: %.1e\n", t_HSS_construct);

        hss_rank = G_HSS.Rank();
        fprintf("  HSS rank: %d\n", hss_rank);
        hss_mem = G_HSS.Storage();
        ratio = hss_mem / N^2;
        fprintf("  ratio: %.1e\n", ratio);

        tic;
        G_HSS.ULV_Factor();
        t_HSS_factor = toc;
        fprintf("  t_HSS_factor: %.1e\n", t_HSS_factor);

        % Direct solution.
        fprintf("Direct solution. \n");
        tic;
        f_direct = G_HSS.ULV_Solve(apply_bf_adj(K_BF, Kf));
        t_solve_direct = toc;

        fprintf("  t_solve_direct: %.1e\n", t_solve_direct);
        rel_res_direct = norm(Kf - apply_bf(K_BF, f_direct)) / norm(Kf);
        fprintf("  rel_res_direct: %.1e\n", rel_res_direct);
        rel_err_direct = norm(f_ex - f_direct) / norm(f_ex);
        fprintf("  rel_res_direct: %.1e\n", rel_err_direct);

        % Iterative solution.
        fprintf("iterative solution.\n");
        rhs = apply_bf_adj(K_BF, Kf);

        fprintf("  without precond.\n");
        tic;
        [f_cg, flag_cg, ~, iter_cg] = pcg(op_G, rhs, tol_cg, maxit_cg);
        t_cg = toc;

        fprintf("    t_cg: %.1e\n", t_cg);
        fprintf("    iter_cg: %d\n", iter_cg);
        rel_res_cg = norm(Kf - apply_bf(K_BF, f_cg)) / norm(Kf);
        fprintf("    rel_res_cg: %.1e\n", rel_res_cg);
        rel_err_cg = norm(f_ex - f_cg) / norm(f_ex);
        fprintf("    rel_err_cg: %.1e\n", rel_err_cg);

        fprintf("  with precond.\n");
        op_M = @(v) G_HSS.ULV_Solve(v);
        tic;
        [f_pcg, flag_pcg, ~, iter_pcg] = pcg(op_G, rhs, tol_cg, maxit_cg, op_M);
        t_pcg = toc;

        fprintf("    t_pcg: %.1e\n", t_pcg);
        fprintf("    iter_pcg: %d\n", iter_pcg);
        rel_res_pcg = norm(Kf - apply_bf(K_BF, f_pcg)) / norm(Kf);
        fprintf("    rel_res_pcg: %.1e\n", rel_res_pcg);
        rel_err_pcg = norm(f_ex - f_pcg) / norm(f_ex);
        fprintf("    rel_err_pcg: %.1e\n", rel_err_pcg);

        save("./data/1d_const_amplituide_results_" + string(p) + "_" + string(tol_hss) + ".mat", ...
            "N", ...
            ...
            "t_BF_construct", ...
            ...
            "t_BF_apply", "rel_err_BF", ...
            ...
            "t_HSS_construct", "hss_rank", "hss_mem", ...
            ...
            "t_HSS_factor", ...
            ...
            "t_solve_direct", "rel_res_direct", "rel_err_direct", ...
            ...
            "t_cg", "iter_cg", "rel_res_cg", "rel_err_cg", ...
            ...
            "t_pcg", "iter_pcg", "rel_res_pcg", "rel_err_pcg");
        clear K_BF;
        clear G_HSS;
        clear f_ex;
        clear Kf;
        clear f_direct;
        clear f_cg;
        clear f_pcg;
    end
end