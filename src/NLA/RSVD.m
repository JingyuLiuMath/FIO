function [U, S, V] = RSVD(A, rank_or_tol, power_iter)
% RSVD Randomized SVD.

% Jingyu Liu, January 31, 2024.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
    power_iter (1, 1) double = 0;
end

arguments (Output)
    U (:, :) double;
    S (:, :) double;
    V (:, :) double;
end

if isempty(A)
    % m = 0; n = 0;
    U = zeros(0, 0);
    S = zeros(0, 0);
    V = zeros(0, 0);
    return;
end

[Q, B] = RangeFinder(A, rank_or_tol, power_iter);
[U1, S, V] = svd(B, "econ");
if rank_or_tol >= 1
    k = min(rank_or_tol, min(size(A)));
    U = Q * U1(:, 1 : k);
    S = S(1 : k, 1 : k);
    V = V(:, 1 : k);
else
    U = Q * U1;
end

end