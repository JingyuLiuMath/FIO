function [U, S, V] = RSVDImplict(m, n, ...
    AX_fun, XA_fun, rank_or_tol, power_iter)
% RSVD Randomized SVD.

% Jingyu Liu, April 12, 2024.

% If rank_or_tol >= 1, it is treated as target rank. Otherwise it is
% treated as relative tolerance.

arguments (Input)
    m (1, 1) double;
    n (1, 1) double;
    AX_fun function_handle;
    XA_fun function_handle;
    rank_or_tol (1, 1) double;
    power_iter (1, 1) double = 0;
end

arguments (Output)
    U (:, :) double;
    S (:, :) double;
    V (:, :) double;
end

if m == 0 || n == 0
    U = zeros(m, 0);
    S = zeros(0, 0);
    V = zeros(n, 0);
    return;
end

[Q, B] = RangeFinderImplict(m, n, ...
    AX_fun, XA_fun, rank_or_tol, power_iter);
[U1, S, V] = svd(B, "econ");
if rank_or_tol >= 1
    k = min([rank_or_tol, m, n]);
    U = Q * U1(:, 1 : k);
    S = S(1 : k, 1 : k);
    V = V(:, 1 : k);
else
    U = Q * U1;
end

end