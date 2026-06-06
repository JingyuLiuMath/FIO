function L = EvalLagrange(x_inter, x_ev)

arguments (Input)
    x_inter (:, 1)
    x_ev (:, 1);
end

arguments (Output)
    L (:, :) double;
end

n = size(x_ev, 1);
r = size(x_inter, 1);

L = ones(n, r);
for j = 1 : r
    ind = [1 : (j - 1), (j + 1) : r];
    curr = (x_ev - x_inter(ind).') ./ (x_inter(j) - x_inter(ind).');
    L(:, j) = prod(curr, 2);
end


end