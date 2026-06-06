function [y, y1, y2] = SpaceCt(B)

arguments (Input)
    B BFNode2D;
end

arguments (Output)
    y (1, 2) double;
    y1 (1, 1) double;
    y2 (1, 1) double;
end

y1 = B.node1_.SpaceCt();
y2 = B.node2_.SpaceCt();
y = [y1, y2];

end