function eta = FreqCt(B)

arguments (Input)
    B BFNode;
end

arguments (Output)
    eta (1, 1) double;
end

eta = (B.ind_start_ + B.ind_end_) / 2 - B.half_N_;

end