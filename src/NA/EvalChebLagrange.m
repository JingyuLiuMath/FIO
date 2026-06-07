function L = EvalChebLagrange(r, x_ev, a, b, type)

arguments (Input)
    r (1,1) double;
    x_ev (:, 1) double;
    a (1, 1) double = -1;
    b (1, 1) double = 1;
    type (1, 1) double = 2;
end

arguments (Output)
    L (:,:) double
end

if r == 1
    L = ones(length(x_ev),1);
    return
end

mid = (a + b) / 2;
half_length = (b - a) / 2;
x_ev = (x_ev - mid) / half_length;
n = r - 1;

switch type
    case 1
        j = (n:-1:0)';
        theta = (2*j + 1) * pi / (2*r);
        x_interp = cos(theta);
        w = (-1).^j .* sin(theta);
    case 2
        j = (n:-1:0)';
        theta = j * pi / n;
        x_interp = cos(theta);
        w = (-1).^j;
        w([1, r]) = w([1, r]) / 2;
end

D = x_ev - x_interp.';

tol = 1e-14;
mask = abs(D) < tol;

A = w.' ./ D;
A(mask) = 0;
L = A ./ sum(A,2);

% Correction
[row, col] = find(mask);
if ~isempty(row)
    rows = unique(row);
    L(rows, :) = 0;
    idx = sub2ind(size(L), row, col);
    L(idx) = 1;
end

end