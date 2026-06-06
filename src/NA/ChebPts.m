function pts = ChebPts(r, a, b)

arguments (Input)
    r (1, 1) double;
    a (1, 1) double = -1;
    b (1, 1) double = 1;
end

arguments (Output)
    pts (:, 1) double;
end

tt = linspace(pi, 0, r)';
pts = cos(tt);
pts = ((b - a) * pts + (b + a)) / 2;

end