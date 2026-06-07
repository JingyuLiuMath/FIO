function L = SpaceEvalChebLagrange(B, r, x_ev)

arguments (Input)
    B BFNode;
    r (1, 1) double;
    x_ev (:, 1) double;
end

arguments (Output)
    L (:, :) double;
end

hh = 1 / B.N_ / 2;
L = EvalChebLagrange(...
    r, ...
    x_ev, ...
    B.ind_start_ / B.N_ - hh, ...
    B.ind_end_ / B.N_ + hh);

end