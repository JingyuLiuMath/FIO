function [kappa1_est, info] = Condest_Implict(op_A, op_A_inv, n, maxit, restarts)
% CONDEST_IMPLICT Estimate 1-norm condition number with Hager algorithm.
%
% Usage:
%   c = Condest_Implict(op_A, op_A_inv, n)
%   [c, info] = Condest_Implict(op_A, op_A_inv, n, maxit, restarts)
%
% Required operators:
%   op_A      : apply A and A^H to vectors.
%   op_A_inv  : apply A^{-1} and A^{-H} to vectors.
%
% Supported call signature for each operator:
%   y = op(x, 'notransp') and y = op(x, 'transp')
%
% Output:
%   kappa1_est : estimate of cond_1(A) = ||A||_1 * ||A^{-1}||_1
%   info       : diagnostics for both Hager iterations.

% Jingyu Liu, April 5, 2026.

arguments (Input)
    op_A (1, 1) function_handle;
    op_A_inv (1, 1) function_handle;
    n (1, 1) double;
    maxit (1, 1) double {mustBeInteger, mustBePositive} = 8;
    restarts (1, 1) double {mustBeInteger, mustBePositive} = 2;
end

arguments (Output)
    kappa1_est (1, 1) double;
    info (1, 1) struct;
end

[normA1_est, infoA] = hager_normest1(op_A, n, maxit, restarts);
[normInvA1_est, infoInvA] = hager_normest1(op_A_inv, n, maxit, restarts);

kappa1_est = normA1_est * normInvA1_est;

info = struct();
info.n = n;
info.normA1_est = normA1_est;
info.normInvA1_est = normInvA1_est;
info.kappa1_est = kappa1_est;
info.A = infoA;
info.A_inv = infoInvA;

end

function [est, info] = hager_normest1(op, n, maxit, restarts)
% Hager's 1-norm estimator using forward and transpose operator actions.

est = 0;
best_x = [];
iter_total = 0;
terminated_by_repeat = false;

for r = 1 : restarts
    if r == 1
        x = ones(n, 1) / n;
    else
        x = randn(n, 1);
        x = x / norm(x, 1);
    end

    used_idx = zeros(maxit, 1);
    used_count = 0;

    for k = 1 : maxit
        iter_total = iter_total + 1;

        y = op(x, 'notransp');
        est_k = norm(y, 1);
        if est_k > est
            est = est_k;
            best_x = x;
        end

        s = sign(y);
        s(s == 0) = 1;
        z = op(s, 'transp');

        [zmax, j] = max(abs(z));
        hx = real(z' * x);

        % Hager stopping criterion: no coordinate gives further ascent.
        if zmax <= hx
            break;
        end

        if any(used_idx(1 : used_count) == j)
            terminated_by_repeat = true;
            break;
        end

        used_count = used_count + 1;
        used_idx(used_count) = j;

        x = zeros(n, 1);
        x(j) = 1;
    end
end

if isempty(best_x)
    est = 0;
end

info = struct();
info.iter_total = iter_total;
info.maxit = maxit;
info.restarts = restarts;
info.terminated_by_repeat = terminated_by_repeat;
info.best_x = best_x;

end