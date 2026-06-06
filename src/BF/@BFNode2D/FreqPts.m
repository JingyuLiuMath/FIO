function [pts, pts1, pts2] = FreqPts(B)

arguments (Input)
    B BFNode2D;
end

arguments (Output)
    pts (:, 2) double;
    pts1 (:, 1) double;
    pts2 (:, 1) double;
end

pts1 = B.node1_.FreqPts();
pts2 = B.node2_.FreqPts();
pts = TensorProduct2D(pts1, pts2);

end