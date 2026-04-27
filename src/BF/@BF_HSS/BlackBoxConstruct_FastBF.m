function BlackBoxConstruct_FastBF(G, K_BF, r, tol)

arguments (Input)
    G BF_HSS;
    K_BF;
    r (1, 1) double;
    tol (1, 1) double;
end

op_G = @(v) apply_fbf_adj(K_BF, apply_fbf(K_BF, v));
G.BlackBoxConstruct(op_G, r, tol);

end