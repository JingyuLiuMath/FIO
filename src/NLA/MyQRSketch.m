function [Q, R, p, k] = MyQRSketch(A, rank_or_tol)
% MyQRSketch For a m-by-n matrix A, compute A(:, p) = Q * R where Q is
% m-by-k and R is k-by-n.

% Jingyu Liu, January 5, 2024.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
end

arguments (Output)
    Q (:, :) double;
    R (:, :) double;
    p (1, :) double;
    k (1, 1) double;
end

if isempty(A)
    % m = 0; n = 0;
    Q = zeros(0, 0);
    R = zeros(0, 0);
    p = [];
    k = 0;
    return;
end

% A(:, p) = Q * R where abs(diag(R)) is is decreasing.
[Q, R, p] = qr(A, "econ", "vector");

if rank_or_tol >= 1
    k = min(rank_or_tol, size(Q, 2));
    Q = Q(:, 1 : k);
    R = R(1 : k, :);
else
    k1 = find(abs(diag(R)) >= rank_or_tol * 1e-1 * max(abs(R(1, 1)), 1), ...
        1, "last");
    k = min(size(R, 1));
    if ~isempty(k1)
        k = min(k1, k);
    else
        k = 0;
    end
    Q = Q(:, 1 : k);
    R = R(1 : k, :);
end

if nargout < 3
    [~, p_inv] = sort(p);
    R = R(:, p_inv);
end

end