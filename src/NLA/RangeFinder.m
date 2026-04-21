function [Q, B] = RangeFinder(A, rank_or_tol, power_iter)
% RangeFinder
% rank_or_tol >= 1:
% Find orthogonal Q (m-by-l) and B (l-by-n) of a m-by-n A such that Q
% approximates the leading l left singular vectors. 
% l = k + p where k = rank_or_tol and p = oversampling_number.
% rank_or_tol < 1:
% Find orthonormal Q (m-by-l) and B (l-by-n) of a m-by-n matrix A such that
% norm(A - Q * B, "fro") / norm(A, "fro") < rank_or_tol.

% Jingyu Liu, January 6, 2024.

arguments (Input)
    A (:, :) double;
    rank_or_tol (1, 1) double;
    power_iter (1, 1) double = 0;
end

arguments (Output)
    Q (:, :) double;
    B (:, :) double;
end

[m, n] = size(A);
min_szA = min(m, n);

if rank_or_tol >= 1
    oversampling_number = min(min_szA, 5);
    k = min(rank_or_tol, min_szA);
    Omega = randn(n, k + oversampling_number);
    [Q, ~] = qr(A * Omega, "econ");
    for it = 1 : power_iter
        [Q, ~] = qr(A' * Q, "econ");
        [Q, ~] = qr(A * Q, "econ");
    end
    % Q = Q(:, 1 : k);
    B = Q' * A;
else
    % rank_or_tol = max(rank_or_tol, sqrt(eps));
    
    max_subspace_dim = max(min_szA, 1);
    block_size = min(max(floor(0.1 * min_szA), 5), min_szA);
    max_iter = ceil(max_subspace_dim / block_size);
    
    block_size0 = block_size;
    Q = zeros(m, 0);
    B = zeros(0, n);
    normA = norm(A, "fro");
    normB = 0;
    old_approx_err = Inf;

    % All qr_econ can be replaced by orth.
    for i = 1 : max_iter
        Omega_i = randn(n, block_size);
        [Q_i, ~] = qr(A * Omega_i - Q * (B * Omega_i), "econ");

        % Power iterations to increase accuracy
        for it = 1 : power_iter
            [Q_i, ~] = qr(A' * Q_i - B' *(Q' * Q_i), "econ");
            [Q_i, ~] = qr(A * Q_i - Q * (B * Q_i), "econ");
        end

        % Reorthogonalization.
        [Q_i, ~] = qr(Q_i - Q * (Q' * Q_i), "econ");

        B_i = Q_i' * A;
        Q = [Q, Q_i];
        B = [B; B_i];
        B_dim = size(B, 1);

        normB = hypot(normB, norm(B_i, "fro"));
        % approx_err(i) = sqrt(abs(normA - normB) * (normA + normB)) / normA;
        approx_err = realsqrt(abs(normA - normB) * (normA + normB)) / normA;

        if approx_err < rank_or_tol
            break;
        end

        if approx_err > old_approx_err
            % Ignore the last step.
            Q(:, (end - block_size + 1) : end) = [];
            B((end - block_size + 1) : end, :) = [];
            break
        end

        if approx_err > old_approx_err / 2
            % Increase block_size.
            block_size = min(block_size + block_size0, ...
                max(max_subspace_dim - B_dim, 1));
        end

        % Make sure max_subspace_dim is not exceeded on last iteration
        if B_dim + block_size > max_subspace_dim
            block_size = max_subspace_dim - B_dim;
        end

        % If bSize is set to 0, break
        if block_size == 0
            break
        end

        old_approx_err = approx_err;
    end
end

end