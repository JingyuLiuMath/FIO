function [L, L1, L2] = SpaceEvalChebLagrange(B, r, x1_ev, x2_ev)

arguments (Input)
    B BFNode2D;
    r (1, 1) double;
    x1_ev (:, 1) double;
    x2_ev (:, 1) double;
end

arguments (Output)
    L (:, :) double;
    L1 (:, :) double;
    L2 (:, :) double;
end

L1 = B.node1_.SpaceEvalChebLagrange(r, x1_ev);
L2 = B.node2_.SpaceEvalChebLagrange(r, x2_ev);
L = kron(L2, L1);

end