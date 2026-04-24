function f = MyApply(G, rhs)

arguments (Input)
    G BF_HSS2D;
    rhs (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

rhs = rhs(G.perm_, :);
f = G.Apply(rhs);
f = f(G.perm_inv_, :);

end