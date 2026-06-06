function [pts, pts1, pts2] = SpacePts(B)

arguments (Input)
    B BFNode2D;
end

arguments (Output)
    pts (:, 2) double;
    pts1 (:, 1) double;
    pts2 (:, 1) double;
end

pts1 = B.node1_.SpacePts();
pts2 = B.node2_.SpacePts();
pts = TensorProduct2D(pts1, pts2);

end