function [U, S, V] = MySVDSketch(A, rank_or_tol)
% MySVDSketch

% Jingyu Liu, May 9, 2024.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    U (:, :) double;
    S (:, :) double;
    V (:, :) double;
end

if rank_or_tol >= 1
    [U, S, V] = svds(A, rank_or_tol);
else
    [U, S, V] = svd(A, "econ");
    k = find(diag(S) >= rank_or_tol * S(1, 1), 1, "last");
    if isempty(k)
        k = 0;
    end
    U = U(:, 1 : k);
    S = S(1 : k, 1 : k);
    V = V(:, 1 : k);
end

end