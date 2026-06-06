function [pts, pts1, pts2] = SpaceChebPts(B, r)

arguments (Input)
    B BFNode2D;
    r (1, 1) double;
end

arguments (Output)
    pts (:, 2) double;
    pts1 (:, 1) double;
    pts2 (:, 1) double;
end

pts1 = B.node1_.SpaceChebPts(r);
pts2 = B.node2_.SpaceChebPts(r);
pts = TensorProduct2D(pts1, pts2);

end