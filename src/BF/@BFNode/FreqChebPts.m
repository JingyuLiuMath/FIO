function pts = FreqChebPts(B, r)

arguments (Input)
    B BFNode;
    r (1, 1) double;
end

arguments (Output)
    pts (:, 1) double;
end

hh = 1 / 2;
pts = ChebPts(r, B.ind_start_ - hh - B.half_N_, B.ind_end_ + hh - B.half_N_);

end