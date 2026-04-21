function [Q, k] = ColBasis(B, k, tol)
% ColBasis

% Jingyu Liu, December 4, 2024.

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
    % m = 0; n = 0;
    Q = zeros(0, 0);
    k = 0;
    return;
end

[Q, R, ~] = qr(B, "econ", "vector");

k1 = find(abs(diag(R)) >= tol * 1e-1 * max(abs(R(1, 1)), 1), 1, "last");
if ~isempty(k1)
    k = min(k1, k);
else
    k = 0;
end

Q = Q(:, 1 : k); 

end