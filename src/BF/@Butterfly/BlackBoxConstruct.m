function out = BlackBoxConstruct(...
    BF, op_K, op_K_adjoint, target_rank, tol, verbose)
% BlackBoxConstruct constructs a BHP factorization from matrix-vector products.

arguments (Input)
    BF Butterfly;
    op_K function_handle;
    op_K_adjoint function_handle;
    target_rank (1, 1) double;
    tol (1, 1) double;
    verbose (1, 1) double = 1;
end

arguments (Output)
    out (1, 1) struct;
end

t_total_start = tic;

target_rank = ceil(target_rank);
BF.constructed_ = false;
assert(target_rank >= 1, "Butterfly:BlackBoxConstruct:InvalidRank", ...
    "target_rank must be positive.");
assert(tol >= 0, "Butterfly:BlackBoxConstruct:InvalidTolerance", ...
    "tol must be nonnegative.");

p = 5;
num_samples = min(BF.N_, target_rank + p);
max_batch_columns = 256;

if verbose == 1
    fprintf("  Black-box BHP construction.\n");
    fprintf("    target rank: %d\n", target_rank);
    fprintf("    num samples per query: %d\n", num_samples);
end

% Construct independent bases for all complementary levels.  A bounded
% batch of block-column sketches shares each operator call.
[left_basis, num_apply_K, middle_left_sketch, middle_omega, ...
    stats_left] = ...
    ConstructLeftBases(...
    BF, op_K, target_rank, tol, num_samples, ...
    max_batch_columns, verbose);
[right_basis, num_apply_K_adjoint, stats_right] = ...
    ConstructRightBases(...
    BF, op_K_adjoint, target_rank, tol, num_samples, ...
    max_batch_columns, verbose);

% Convert the independently sampled bases directly to BHP transfer blocks.
t_left_factor_start = tic;
ConstructLeftFactors(BF, left_basis);
out.left_factor_time = toc(t_left_factor_start);
t_right_factor_start = tic;
ConstructRightFactors(BF, right_basis);
out.right_factor_time = toc(t_right_factor_start);

% Recover the middle coupling matrices from the projected sketches saved
% during the middle-level left-basis construction.
t_middle_start = tic;
stats_middle = ConstructMiddleFactors(...
    BF, middle_left_sketch, middle_omega, right_basis{1}, tol);
out.middle_time = toc(t_middle_start);

BF.constructed_ = true;

out.left_total_time = stats_left.total_time;
out.left_operator_time = stats_left.operator_time;
out.left_basis_time = stats_left.basis_time;
out.left_projection_time = stats_left.projection_time;
out.right_total_time = stats_right.total_time;
out.right_operator_time = stats_right.operator_time;
out.right_basis_time = stats_right.basis_time;
out.middle_pinv_time = stats_middle.pinv_time;
out.operator_time = ...
    out.left_operator_time + out.right_operator_time;
out.basis_time = out.left_basis_time + out.right_basis_time;
out.factor_time = out.left_factor_time + out.right_factor_time;
out.num_apply_K = num_apply_K;
out.num_apply_K_adjoint = num_apply_K_adjoint;
out.num_operator_calls_K = stats_left.num_operator_calls;
out.num_operator_calls_K_adjoint = stats_right.num_operator_calls;
out.max_batch_columns = max_batch_columns;
out.total_time = toc(t_total_start);
out.other_time = max(0, out.total_time ...
    - out.operator_time - out.basis_time ...
    - out.left_projection_time - out.factor_time ...
    - out.middle_time);

if verbose == 1
    fprintf("    num K applications: %d\n", num_apply_K);
    fprintf("    num K' applications: %d\n", num_apply_K_adjoint);
    fprintf("    num K calls: %d\n", out.num_operator_calls_K);
    fprintf("    num K' calls: %d\n", ...
        out.num_operator_calls_K_adjoint);
end

end

function [basis_level, num_apply, middle_left_sketch, middle_omega, ...
    stats] = ...
    ConstructLeftBases(...
    BF, op_K, target_rank, tol, num_samples, ...
    max_batch_columns, verbose)

t_total_start = tic;
stats.operator_time = 0;
stats.basis_time = 0;
stats.projection_time = 0;
stats.num_operator_calls = 0;

L = BF.L_;
h_x = BF.h_x_;
L_x = BF.L_x_;
N = BF.N_;

basis_level = cell(1, L_x - h_x + 1);
num_apply = 0;
middle_row_nodes = BF.tree_{h_x + 1};
middle_col_nodes = BF.tree_{L - h_x + 1};
middle_left_sketch = cell(...
    length(middle_row_nodes), length(middle_col_nodes));
middle_omega = cell(1, length(middle_col_nodes));
max_blocks_per_batch = max(1, floor(max_batch_columns / num_samples));
Omega = zeros(N, max_blocks_per_batch * num_samples);

for level_x = h_x : L_x
    level_xi = L - level_x;
    row_nodes = BF.tree_{level_x + 1};
    col_nodes = BF.tree_{level_xi + 1};
    m_x = length(row_nodes);
    m_xi = length(col_nodes);

    basis = cell(m_x, m_xi);
    for batch_start = 1 : max_blocks_per_batch : m_xi
        batch_end = min(batch_start + max_blocks_per_batch - 1, m_xi);
        num_blocks = batch_end - batch_start + 1;
        num_rhs = num_blocks * num_samples;
        Omega(:, 1 : num_rhs) = 0;

        for ind_sigma = batch_start : batch_end
            sigma = col_nodes{ind_sigma};
            ind_col = (sigma.ind_start_ + 1) : (sigma.ind_end_ + 1);
            ind_batch = ind_sigma - batch_start;
            ind_rhs = ind_batch * num_samples + (1 : num_samples);
            Omega(ind_col, ind_rhs) = ...
                randn(length(ind_col), num_samples);
        end

        t_apply_start = tic;
        Y = ApplyOperator(op_K, Omega(:, 1 : num_rhs), N);
        stats.operator_time = stats.operator_time + toc(t_apply_start);
        stats.num_operator_calls = stats.num_operator_calls + 1;
        num_apply = num_apply + num_rhs;

        for ind_sigma = batch_start : batch_end
            sigma = col_nodes{ind_sigma};
            ind_col = (sigma.ind_start_ + 1) : (sigma.ind_end_ + 1);
            ind_batch = ind_sigma - batch_start;
            ind_rhs = ind_batch * num_samples + (1 : num_samples);
            for ind_tau = 1 : m_x
                tau = row_nodes{ind_tau};
                ind_row = (tau.ind_start_ + 1) : (tau.ind_end_ + 1);
                Y_tau_sigma = Y(ind_row, ind_rhs);
                t_basis_start = tic;
                basis{ind_tau, ind_sigma} = ColBasis(...
                    Y_tau_sigma, target_rank, tol);
                stats.basis_time = stats.basis_time + toc(t_basis_start);
                if level_x == h_x
                    t_projection_start = tic;
                    middle_left_sketch{ind_tau, ind_sigma} = ...
                        basis{ind_tau, ind_sigma}' * Y_tau_sigma;
                    stats.projection_time = stats.projection_time ...
                        + toc(t_projection_start);
                end
            end
            if level_x == h_x
                middle_omega{ind_sigma} = Omega(ind_col, ind_rhs);
            end
        end
    end
    basis_level{level_x - h_x + 1} = basis;

    if verbose == 1
        fprintf("    left level %d: max rank %d\n", ...
            level_x, MaxCellRank(basis));
    end
end

stats.total_time = toc(t_total_start);

end

function [basis_level, num_apply, stats] = ConstructRightBases(...
    BF, op_K_adjoint, target_rank, tol, num_samples, ...
    max_batch_columns, verbose)

t_total_start = tic;
stats.operator_time = 0;
stats.basis_time = 0;
stats.num_operator_calls = 0;

L = BF.L_;
h_xi = BF.h_xi_;
L_xi = BF.L_xi_;
N = BF.N_;

basis_level = cell(1, L_xi - h_xi + 1);
num_apply = 0;
max_blocks_per_batch = max(1, floor(max_batch_columns / num_samples));
Psi = zeros(N, max_blocks_per_batch * num_samples);

for level_xi = h_xi : L_xi
    level_x = L - level_xi;
    row_nodes = BF.tree_{level_x + 1};
    col_nodes = BF.tree_{level_xi + 1};
    m_x = length(row_nodes);
    m_xi = length(col_nodes);

    basis = cell(m_x, m_xi);
    for batch_start = 1 : max_blocks_per_batch : m_x
        batch_end = min(batch_start + max_blocks_per_batch - 1, m_x);
        num_blocks = batch_end - batch_start + 1;
        num_rhs = num_blocks * num_samples;
        Psi(:, 1 : num_rhs) = 0;

        for ind_tau = batch_start : batch_end
            tau = row_nodes{ind_tau};
            ind_row = (tau.ind_start_ + 1) : (tau.ind_end_ + 1);
            ind_batch = ind_tau - batch_start;
            ind_rhs = ind_batch * num_samples + (1 : num_samples);
            Psi(ind_row, ind_rhs) = ...
                randn(length(ind_row), num_samples);
        end

        t_apply_start = tic;
        Z = ApplyOperator(op_K_adjoint, Psi(:, 1 : num_rhs), N);
        stats.operator_time = stats.operator_time + toc(t_apply_start);
        stats.num_operator_calls = stats.num_operator_calls + 1;
        num_apply = num_apply + num_rhs;

        for ind_tau = batch_start : batch_end
            ind_batch = ind_tau - batch_start;
            ind_rhs = ind_batch * num_samples + (1 : num_samples);
            for ind_sigma = 1 : m_xi
                sigma = col_nodes{ind_sigma};
                ind_col = (sigma.ind_start_ + 1) : (sigma.ind_end_ + 1);
                Z_tau_sigma = Z(ind_col, ind_rhs);
                t_basis_start = tic;
                basis{ind_tau, ind_sigma} = ColBasis(...
                    Z_tau_sigma, target_rank, tol);
                stats.basis_time = stats.basis_time + toc(t_basis_start);
            end
        end
    end
    basis_level{level_xi - h_xi + 1} = basis;

    if verbose == 1
        fprintf("    right level %d: max rank %d\n", ...
            level_xi, MaxCellRank(basis));
    end
end

stats.total_time = toc(t_total_start);

end

function ConstructLeftFactors(BF, basis_level)

L = BF.L_;
h_x = BF.h_x_;
L_x = BF.L_x_;
num_children = BF.num_children_;
ch_list = BF.ch_list_;

BF.G_ = cell(1, L_x - h_x);
for level_x = h_x : (L_x - 1)
    level_xi = L - level_x;
    prev_basis = basis_level{level_x - h_x + 1};
    curr_basis = basis_level{level_x - h_x + 2};

    row_nodes = BF.tree_{level_x + 1};
    row_child_nodes = BF.tree_{level_x + 2};
    col_parent_nodes = BF.tree_{level_xi};

    m_x = length(row_child_nodes);
    m_xi = length(col_parent_nodes);
    cnt_G = L_x - level_x;
    BF.G_{cnt_G} = cell(m_x, m_xi);

    for ind_alpha = 1 : length(row_nodes)
        alpha = row_nodes{ind_alpha};
        for ind_sigma = 1 : m_xi
            sigma = col_parent_nodes{ind_sigma};
            for ch_alpha = ch_list
                ind_tau = num_children * alpha.order_ + ch_alpha + 1;
                tau = row_child_nodes{ind_tau};
                ind_tau_local = ...
                    (tau.ind_start_ - alpha.ind_start_ + 1) : ...
                    (tau.ind_end_ - alpha.ind_start_ + 1);

                U_curr = curr_basis{ind_tau, ind_sigma};
                beta_rank = zeros(1, num_children);
                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    U_prev = prev_basis{ind_alpha, ind_beta};
                    beta_rank(ch_sigma + 1) = size(U_prev, 2);
                end

                G_tau_sigma = zeros(size(U_curr, 2), ...
                    sum(beta_rank), "like", U_curr);
                beta_offset = 0;
                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    U_prev = prev_basis{ind_alpha, ind_beta};
                    rank_beta = beta_rank(ch_sigma + 1);
                    G_tau_sigma(:, ...
                        (beta_offset + 1) : (beta_offset + rank_beta)) = ...
                        U_curr' * U_prev(ind_tau_local, :);
                    beta_offset = beta_offset + rank_beta;
                end
                BF.G_{cnt_G}{ind_tau, ind_sigma} = G_tau_sigma;
            end
        end
    end
end

BF.U_ = basis_level{end};

end

function ConstructRightFactors(BF, basis_level)

L = BF.L_;
h_xi = BF.h_xi_;
L_xi = BF.L_xi_;
num_children = BF.num_children_;
ch_list = BF.ch_list_;

BF.H_ = cell(1, L_xi - h_xi);
for level_xi = h_xi : (L_xi - 1)
    level_x = L - level_xi;
    prev_basis = basis_level{level_xi - h_xi + 1};
    curr_basis = basis_level{level_xi - h_xi + 2};

    row_parent_nodes = BF.tree_{level_x};
    col_nodes = BF.tree_{level_xi + 1};
    col_child_nodes = BF.tree_{level_xi + 2};

    m_x = length(row_parent_nodes);
    m_xi = length(col_child_nodes);
    cnt_H = level_xi - h_xi + 1;
    BF.H_{cnt_H} = cell(m_x, m_xi);

    for ind_beta = 1 : length(col_nodes)
        beta = col_nodes{ind_beta};
        for ind_tau = 1 : m_x
            tau = row_parent_nodes{ind_tau};
            for ch_beta = ch_list
                ind_sigma = num_children * beta.order_ + ch_beta + 1;
                sigma = col_child_nodes{ind_sigma};
                ind_sigma_local = ...
                    (sigma.ind_start_ - beta.ind_start_ + 1) : ...
                    (sigma.ind_end_ - beta.ind_start_ + 1);

                V_curr = curr_basis{ind_tau, ind_sigma};
                alpha_rank = zeros(num_children, 1);
                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    V_prev = prev_basis{ind_alpha, ind_beta};
                    alpha_rank(ch_tau + 1) = size(V_prev, 2);
                end

                H_tau_sigma = zeros(sum(alpha_rank), ...
                    size(V_curr, 2), "like", V_curr);
                alpha_offset = 0;
                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    V_prev = prev_basis{ind_alpha, ind_beta};
                    rank_alpha = alpha_rank(ch_tau + 1);
                    H_tau_sigma(...
                        (alpha_offset + 1) : (alpha_offset + rank_alpha), :) = ...
                        V_prev(ind_sigma_local, :)' * V_curr;
                    alpha_offset = alpha_offset + rank_alpha;
                end
                BF.H_{cnt_H}{ind_tau, ind_sigma} = H_tau_sigma;
            end
        end
    end
end

outer_basis = basis_level{end};
BF.V_ = cell(size(outer_basis));
for i = 1 : size(outer_basis, 1)
    for j = 1 : size(outer_basis, 2)
        BF.V_{i, j} = outer_basis{i, j}';
    end
end

end

function stats = ConstructMiddleFactors(...
    BF, middle_left_sketch, middle_omega, right_basis, tol)

stats.pinv_time = 0;

h_x = BF.h_x_;
h_xi = BF.h_xi_;

row_nodes = BF.tree_{h_x + 1};
col_nodes = BF.tree_{h_xi + 1};
m_x = length(row_nodes);
m_xi = length(col_nodes);

BF.M_ = cell(m_x, m_xi);
for ind_sigma = 1 : m_xi
    Omega_sigma = middle_omega{ind_sigma};

    for ind_tau = 1 : m_x
        V = right_basis{ind_tau, ind_sigma};
        right_sketch = V' * Omega_sigma;
        t_pinv_start = tic;
        right_pinv = StablePinv(right_sketch, tol);
        stats.pinv_time = stats.pinv_time + toc(t_pinv_start);
        BF.M_{ind_tau, ind_sigma} = ...
            middle_left_sketch{ind_tau, ind_sigma} ...
            * right_pinv;
    end
end

end

function A_pinv = StablePinv(A, tol)

if isempty(A)
    A_pinv = zeros(size(A, 2), size(A, 1), "like", A);
    return;
end

[U, S, V] = svd(A, "econ");
s = diag(S);
if isempty(s) || s(1) == 0
    A_pinv = zeros(size(A, 2), size(A, 1), "like", A);
    return;
end

threshold = tol * 1e-1 * s(1);
ind = find(s >= threshold & s > 0);
A_pinv = V(:, ind) * diag(1 ./ s(ind)) * U(:, ind)';

end

function Y = ApplyOperator(op_A, X, N)

try
    Y = op_A(X);
    if ~isequal(size(Y), [N, size(X, 2)])
        error("Butterfly:BlackBoxConstruct:InvalidOperatorOutput", ...
            "The black-box operator returned an array of invalid size.");
    end
catch exception
    if size(X, 2) == 1
        rethrow(exception);
    end

    Y = zeros(N, size(X, 2), "like", X);
    for j = 1 : size(X, 2)
        y = op_A(X(:, j));
        if ~isequal(size(y), [N, 1])
            rethrow(exception);
        end
        Y(:, j) = y;
    end
end

end

function r = MaxCellRank(C)

r = 0;
for i = 1 : numel(C)
    r = max(r, size(C{i}, 2));
end

end
