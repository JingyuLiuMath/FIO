function BlackBoxConstruct_Indep_FastBF(G, K_BF, rank_func, tol)

arguments (Input)
    G BF_HSS2D;
    K_BF;
    rank_func function_handle;
    tol (1, 1) double;
end

op_G = @(v) G_fun(v, K_BF, G.perm_, G.perm_inv_);
G.BlackBoxConstruct_Indep(op_G, rank_func, tol);

end

function y = G_fun(v, K_BF, p, p_inv)

v = v(p_inv, :);
y = apply_fbf_adj(K_BF, apply_fbf(K_BF, v));
y = y(p, :);

end