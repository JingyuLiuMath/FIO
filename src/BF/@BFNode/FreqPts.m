function pts = FreqPts(B)

arguments (Input)
    B BFNode;
end

arguments (Output)
    pts (:, 1) double;
end

pts = (B.ind_start_ : B.ind_end_)' - B.half_N_;

end