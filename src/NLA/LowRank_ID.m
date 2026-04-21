function [sk, V, k, re] = LowRank_ID(A, rank_or_tol)
% LowRank_ID Interpolative decomposition. A = A(:, sk) * V'.

% Jingyu Liu, November 19, 2024.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    sk (1, :) double;
    V (:, :) double;
    k (1, 1) double;
    re (1, :) double;
end

[sk, re, T, k] = ID(A, rank_or_tol);
V = zeros(size(A, 2), length(sk));
V(sk, :) = eye(length(sk));
V(re, :) = T';

% err = norm(A - A(:, sk) * V', "fro") / norm(A, "fro");
% fprintf("err: %.4e\n", err);

end