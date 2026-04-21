function [sk, re, Z, k] = Row_ID(A, rank_or_tol)
% Row_ID Interpolative decomposition. A(re, :) = Z * A(sk, :).

% Jingyu Liu, December 11, 2023.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    sk (1, :) double;
    re (1, :) double;
    Z (:, :) double;
    k (1, 1) double;
end

[sk, re, T, k] = ID(A.', rank_or_tol);
Z = T.';

end