function y = SpaceCt(B)

arguments (Input)
    B BFNode;
end

arguments (Output)
    y (1, 1) double;
end

y = (B.ind_start_ + B.ind_end_) / 2 / B.N_;

end