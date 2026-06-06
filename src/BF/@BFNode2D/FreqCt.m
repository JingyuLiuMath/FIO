function [eta, eta1, eta2] = FreqCt(B)

arguments (Input)
    B BFNode2D;
end

arguments (Output)
    eta (1, 2) double;
    eta1 (1, 1) double;
    eta2 (1, 1) double;
end

eta1 = B.node1_.FreqCt();
eta2 = B.node2_.FreqCt();
eta = [eta1, eta2];

end