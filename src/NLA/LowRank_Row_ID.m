function [sk, U, k, re] = LowRank_Row_ID(A, rank_or_tol)
% LowRank_Row_ID Interpolative decomposition. A = U * A(sk, :).

% Jingyu Liu, November 19, 2024.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    sk (1, :) double;
    U (:, :) double;
    k (1, 1) double;
    re (1, :) double;
end

[sk, re, Z, k] = Row_ID(A, rank_or_tol);
U = zeros(size(A, 1), length(sk));
U(sk, :) = eye(length(sk));
U(re, :) = Z;

% err = norm(A - U * A(sk, :), "fro") / norm(A, "fro");
% fprintf("err: %.4e\n", err);

end