function [L, L1, L2] = FreqEvalChebLagrange(B, r, xi1_ev, xi2_ev)

arguments (Input)
    B BFNode2D;
    r (1, 1) double;
    xi1_ev (:, 1) double;
    xi2_ev (:, 1) double;
end

arguments (Output)
    L (:, :) double;
    L1 (:, :) double;
    L2 (:, :) double;
end

L1 = B.node1_.FreqEvalChebLagrange(r, xi1_ev);
L2 = B.node2_.FreqEvalChebLagrange(r, xi2_ev);
L = kron(L2, L1);

end