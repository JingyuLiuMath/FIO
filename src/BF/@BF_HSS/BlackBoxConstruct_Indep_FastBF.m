function BlackBoxConstruct_Indep_FastBF(G, op_G, rank_func, tol)

arguments (Input)
    G BF_HSS;
    op_G;
    rank_func function_handle;
    tol (1, 1) double;
end

op_G_HSS = @(v) op_G(v);
G.BlackBoxConstruct_Indep(op_G_HSS, rank_func, tol);

end