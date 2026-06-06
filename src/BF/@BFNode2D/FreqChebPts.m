function [pts, pts1, pts2] = FreqChebPts(B, r)

arguments (Input)
    B BFNode2D;
    r (1, 1) double;
end

arguments (Output)
    pts (:, 2) double;
    pts1 (:, 1) double;
    pts2 (:, 1) double;
end


pts1 = B.node1_.FreqChebPts(r);
pts2 = B.node2_.FreqChebPts(r);
pts = TensorProduct2D(pts1, pts2);

end