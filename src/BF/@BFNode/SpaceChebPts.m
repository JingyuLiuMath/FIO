function pts = SpaceChebPts(B, r)

arguments (Input)
    B BFNode;
    r (1, 1) double;
end

arguments (Output)
    pts (:, 1) double;
end

hh = 1 / B.N_ / 2;
pts = ChebPts(r, B.ind_start_ / B.N_ - hh, B.ind_end_ / B.N_ + hh);

end