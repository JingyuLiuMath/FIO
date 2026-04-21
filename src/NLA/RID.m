function [sk, re, T, k] = RID(A, rank_or_tol)
% RID Randomized ID. A(:, re) = A(:, sk) * T

% Jingyu Liu, March 20, 2024.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    sk (1, :) double;
    re (1, :) double;
    T (:, :) double;
    k (1, 1) double;
end

min_szA = min(size(A));
if min_szA <= 1024
    [sk, re, T, k] = ID(A, rank_or_tol);
else
    [~, B] = RangeFinder(A, rank_or_tol);
    [sk, re, T, k] = ID(B, rank_or_tol);
end

end