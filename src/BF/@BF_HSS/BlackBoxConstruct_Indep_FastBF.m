function BlackBoxConstruct_Indep_FastBF(G, K_BF, rank_func, tol)

arguments (Input)
    G BF_HSS;
    K_BF;
    rank_func function_handle;
    tol (1, 1) double;
end

op_G = @(v) apply_fbf_adj(K_BF, apply_fbf(K_BF, v));
G.BlackBoxConstruct_Indep(op_G, rank_func, tol);

end