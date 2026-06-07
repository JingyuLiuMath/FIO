function L = EvalLagrange(x_interp, x_ev)

arguments (Input)
    x_interp (:, 1)
    x_ev (:, 1);
end

arguments (Output)
    L (:, :) double;
end

n = size(x_ev, 1);
r = size(x_interp, 1);

L = ones(n, r);
for j = 1 : r
    ind = [1 : (j - 1), (j + 1) : r];
    curr = (x_ev - x_interp(ind).') ./ (x_interp(j) - x_interp(ind).');
    L(:, j) = prod(curr, 2);
end


end