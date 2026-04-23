function f = Solve(G, p, p_inv, rhs)

arguments (Input)
    G BF_HSS2D;
    p (:, 1) double;
    p_inv (:, 1) double;
    rhs (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

rhs = rhs(p, :);
f = G.ULV_Solve(rhs);
f = f(p_inv, :);

end