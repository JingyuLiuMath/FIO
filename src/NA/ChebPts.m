function pts = ChebPts(r, a, b, type)

arguments (Input)
    r (1, 1) double;
    a (1, 1) double = -1;
    b (1, 1) double = 1;
    type (1, 1) double = 2;
end

arguments (Output)
    pts (:, 1) double;
end

switch type
    case 1
        tt = (2 * (r : -1 : 1) - 1)' * pi / 2 / r;
    case 2
        tt = linspace(pi, 0, r)';
end

pts = cos(tt);
pts = ((b - a) * pts + (b + a)) / 2;

end