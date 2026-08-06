function y = apply_fmbf_adj(Factors, x)

if size(x, 2) == 0
    y = x;
    return;
end

y = zeros(size(x), "like", x);

for i = 1 : (size(Factors, 1) - 1)
    y(Factors{i, 2}, :) = ...
        apply_fbf_adj(Factors{i, 1}, x);
end
y(Factors{end, 2}, :) = Factors{end, 1}' * x;

end
