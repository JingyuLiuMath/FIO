function L = FreqEvalChebLagrange(B, r, xi_ev)

arguments (Input)
    B BFNode;
    r (1, 1) double;
    xi_ev (:, 1) double;
end

arguments (Output)
    L (:, :) double;
end

hh = 1 / 2;
L = EvalChebLagrange(...
    r, ...
    xi_ev, ...
    B.ind_start_ - hh - B.half_N_, ...
    B.ind_end_ + hh - B.half_N_);

end