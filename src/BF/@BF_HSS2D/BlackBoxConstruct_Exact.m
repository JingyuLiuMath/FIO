function BlackBoxConstruct_Exact(G, K, r, tol)

arguments (Input)
    G BF_HSS2D;
    K;
    r (1, 1) double;
    tol (1, 1) double;
end

op_G = @(v) G_fun(v, K, G.perm_, G.perm_inv_);
G.BlackBoxConstruct(op_G, r, tol);

end

function y = G_fun(v, K, p, p_inv)

v = v(p_inv, :);
y = K' * (K * v);
y = y(p, :);

end