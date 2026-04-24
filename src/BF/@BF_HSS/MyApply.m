function f = MyApply(G, rhs)

arguments (Input)
    G BF_HSS;
    rhs (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

f = G.Apply(rhs);

end