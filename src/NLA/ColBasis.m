function [Q, k] = ColBasis(B, k, tol)
% ColBasis

arguments (Input)
    B (:, :) double;
    k (1, 1) double;
    tol (1, 1) double;
end

arguments (Output)
    Q (:, :) double;
    k (1, 1) double;
end

if isempty(B)
    Q = zeros(size(B, 1), 0, "like", B);
    k = 0;
    return;
end

[Q, R, ~] = qr(B, "econ", "vector");

diag_R = abs(diag(R));
if isempty(diag_R) || diag_R(1) == 0
    Q = zeros(size(B, 1), 0, "like", B);
    k = 0;
    return;
end

k1 = find(diag_R >= tol * 1e-1 * diag_R(1), 1, "last");
if ~isempty(k1)
    k = min(k1, k);
else
    k = 0;
end

Q = Q(:, 1 : k); 

end
