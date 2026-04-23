function f = Solve(G, rhs)

arguments (Input)
    G BF_HSS;
    rhs (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

f = G.ULV_Solve(rhs);

end