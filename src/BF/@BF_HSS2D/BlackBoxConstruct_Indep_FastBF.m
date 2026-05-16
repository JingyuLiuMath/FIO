function BlackBoxConstruct_Indep_FastBF(G, op_G, rank_func, tol)

arguments (Input)
    G BF_HSS2D;
    op_G;
    rank_func function_handle;
    tol (1, 1) double;
end

op_G_HSS = @(v) G_fun(v, op_G, G.perm_, G.perm_inv_);
G.BlackBoxConstruct_Indep_New(op_G_HSS, rank_func, tol);

end

function y = G_fun(v, op_G, p, p_inv)

v = v(p_inv, :);
y = op_G(v);
y = y(p, :);

end