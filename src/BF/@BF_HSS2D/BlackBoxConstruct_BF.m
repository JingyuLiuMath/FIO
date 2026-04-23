function BlackBoxConstruct_BF(G, p, p_inv, K_BF, r, tol)

arguments (Input)
    G BF_HSS2D;
    p (:, 1) double;
    p_inv (:, 1) double;
    K_BF;
    r (1, 1) double;
    tol (1, 1) double;
end

op_G = @(v) G_fun(v, K_BF, p, p_inv);
G.BlackBoxConstruct(op_G, r, tol);

end

function y = G_fun(v, K_BF, p, p_inv)

v = v(p_inv, :);
y = apply_fbf_adj(K_BF, apply_fbf(K_BF, v));
y = y(p, :);

end