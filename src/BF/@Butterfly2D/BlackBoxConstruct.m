function out = BlackBoxConstruct( ...
    BF, op_K, op_K_adjoint, target_rank, tol, verbose)
% BlackBoxConstruct constructs a BHP factorization from matrix-vector products.

arguments (Input)
    BF Butterfly2D;
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
assert(target_rank >= 1, "Butterfly2D:BlackBoxConstruct:InvalidRank", ...
    "target_rank must be positive.");
assert(tol >= 0, "Butterfly2D:BlackBoxConstruct:InvalidTolerance", ...
    "tol must be nonnegative.");

p = 5;
num_samples = min(BF.N_, target_rank + p);

if verbose == 1
    fprintf("  Black-box BHP construction.\n");
    fprintf("    target rank: %d\n", target_rank);
    fprintf("    num samples per query: %d\n", num_samples);
end

% Construct independent bases for all complementary levels.  Only one
% block-column sketch is live at a time.
[left_basis, num_apply_K] = ConstructLeftBases(...
    BF, op_K, target_rank, tol, num_samples, verbose);
[right_basis, num_apply_K_adjoint] = ConstructRightBases(...
    BF, op_K_adjoint, target_rank, tol, num_samples, verbose);

% Convert the independently sampled bases directly to BHP transfer blocks.
ConstructLeftFactors(BF, left_basis);
ConstructRightFactors(BF, right_basis);

% Recover the middle coupling matrices.  Recomputing one middle-level
% sketch avoids storing all middle sketches during the basis construction.
num_apply_K = num_apply_K + ConstructMiddleFactors(...
    BF, op_K, left_basis{1}, right_basis{1}, tol, num_samples);

BF.constructed_ = true;
out.num_apply_K = num_apply_K;
out.num_apply_K_adjoint = num_apply_K_adjoint;
out.total_time = toc(t_total_start);

if verbose == 1
    fprintf("    num K applications: %d\n", num_apply_K);
    fprintf("    num K' applications: %d\n", num_apply_K_adjoint);
end

end

function [basis_level, num_apply] = ConstructLeftBases(...
        BF, op_K, target_rank, tol, num_samples, verbose)

L = BF.L_;
h_x = BF.h_x_;
L_x = BF.L_x_;
N = BF.N_;

basis_level = cell(1, L_x - h_x + 1);
num_apply = 0;

for level_x = h_x : L_x
    level_xi = L - level_x;
    row_nodes = BF.tree_{level_x + 1};
    col_nodes = BF.tree_{level_xi + 1};
    m_x = length(row_nodes);
    m_xi = length(col_nodes);

    basis = cell(m_x, m_xi);
    for ind_sigma = 1 : m_xi
        ind_col = NodeIndices(col_nodes{ind_sigma});
        Omega = zeros(N, num_samples);
        Omega(ind_col, :) = GaussianMatrix(length(ind_col), num_samples);
        Y = ApplyOperator(op_K, Omega, N);
        num_apply = num_apply + num_samples;

        for ind_tau = 1 : m_x
            ind_row = NodeIndices(row_nodes{ind_tau});
            basis{ind_tau, ind_sigma} = NumericalBasis(...
                Y(ind_row, :), target_rank, tol);
        end
    end
    basis_level{level_x - h_x + 1} = basis;

    if verbose == 1
        fprintf("    left level %d: max rank %d\n", ...
            level_x, MaxCellRank(basis));
    end
end

end

function [basis_level, num_apply] = ConstructRightBases(...
        BF, op_K_adjoint, target_rank, tol, num_samples, verbose)

L = BF.L_;
h_xi = BF.h_xi_;
L_xi = BF.L_xi_;
N = BF.N_;

basis_level = cell(1, L_xi - h_xi + 1);
num_apply = 0;

for level_xi = h_xi : L_xi
    level_x = L - level_xi;
    row_nodes = BF.tree_{level_x + 1};
    col_nodes = BF.tree_{level_xi + 1};
    m_x = length(row_nodes);
    m_xi = length(col_nodes);

    basis = cell(m_x, m_xi);
    for ind_tau = 1 : m_x
        ind_row = NodeIndices(row_nodes{ind_tau});
        Psi = zeros(N, num_samples);
        Psi(ind_row, :) = GaussianMatrix(length(ind_row), num_samples);
        Z = ApplyOperator(op_K_adjoint, Psi, N);
        num_apply = num_apply + num_samples;

        for ind_sigma = 1 : m_xi
            ind_col = NodeIndices(col_nodes{ind_sigma});
            basis{ind_tau, ind_sigma} = NumericalBasis(...
                Z(ind_col, :), target_rank, tol);
        end
    end
    basis_level{level_xi - h_xi + 1} = basis;

    if verbose == 1
        fprintf("    right level %d: max rank %d\n", ...
            level_xi, MaxCellRank(basis));
    end
end

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
        ind_alpha_global = NodeIndices(alpha);
        for ind_sigma = 1 : m_xi
            sigma = col_parent_nodes{ind_sigma};
            for ch_alpha = ch_list
                ind_tau = num_children * alpha.order_ + ch_alpha + 1;
                tau = row_child_nodes{ind_tau};
                ind_tau_global = NodeIndices(tau);
                ind_tau_local = LocalIndices(ind_alpha_global, ind_tau_global);

                U_curr = curr_basis{ind_tau, ind_sigma};
                G_tau_beta_cell = cell(1, num_children);
                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    U_prev = prev_basis{ind_alpha, ind_beta};
                    G_tau_beta_cell{ch_sigma + 1} = ...
                        U_curr' * U_prev(ind_tau_local, :);
                end
                BF.G_{cnt_G}{ind_tau, ind_sigma} = ...
                    cell2mat(G_tau_beta_cell);
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
        ind_beta_global = NodeIndices(beta);
        for ind_tau = 1 : m_x
            tau = row_parent_nodes{ind_tau};
            for ch_beta = ch_list
                ind_sigma = num_children * beta.order_ + ch_beta + 1;
                sigma = col_child_nodes{ind_sigma};
                ind_sigma_global = NodeIndices(sigma);
                ind_sigma_local = LocalIndices(...
                    ind_beta_global, ind_sigma_global);

                V_curr = curr_basis{ind_tau, ind_sigma};
                H_alpha_sigma_cell = cell(num_children, 1);
                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    V_prev = prev_basis{ind_alpha, ind_beta};
                    H_alpha_sigma_cell{ch_tau + 1} = ...
                        V_prev(ind_sigma_local, :)' * V_curr;
                end
                BF.H_{cnt_H}{ind_tau, ind_sigma} = ...
                    cell2mat(H_alpha_sigma_cell);
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

function num_apply = ConstructMiddleFactors(...
        BF, op_K, left_basis, right_basis, tol, num_samples)

h_x = BF.h_x_;
h_xi = BF.h_xi_;
N = BF.N_;

row_nodes = BF.tree_{h_x + 1};
col_nodes = BF.tree_{h_xi + 1};
m_x = length(row_nodes);
m_xi = length(col_nodes);

BF.M_ = cell(m_x, m_xi);
num_apply = 0;
for ind_sigma = 1 : m_xi
    ind_col = NodeIndices(col_nodes{ind_sigma});
    Omega_sigma = GaussianMatrix(length(ind_col), num_samples);
    Omega = zeros(N, num_samples);
    Omega(ind_col, :) = Omega_sigma;
    Y = ApplyOperator(op_K, Omega, N);
    num_apply = num_apply + num_samples;

    for ind_tau = 1 : m_x
        ind_row = NodeIndices(row_nodes{ind_tau});
        U = left_basis{ind_tau, ind_sigma};
        V = right_basis{ind_tau, ind_sigma};
        right_sketch = V' * Omega_sigma;
        BF.M_{ind_tau, ind_sigma} = ...
            (U' * Y(ind_row, :)) * StablePinv(right_sketch, tol);
    end
end

end

function Q = NumericalBasis(A, target_rank, tol)

[Q, ~] = ColBasis(A, target_rank, tol);

end

function A_pinv = StablePinv(A, tol)

if isempty(A)
    A_pinv = zeros(size(A, 2), size(A, 1));
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
        error("Butterfly2D:BlackBoxConstruct:InvalidOperatorOutput", ...
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

function G = GaussianMatrix(m, n)

G = randn(m, n);

end

function ind = NodeIndices(node)

ind = node.Indices();

end

function ind_local = LocalIndices(ind_parent, ind_child)

[is_member, ind_local] = ismember(ind_child, ind_parent);
assert(all(is_member), "Butterfly2D:BlackBoxConstruct:InvalidTree", ...
    "A child node contains indices outside its parent.");

end

function r = MaxCellRank(C)

r = 0;
for i = 1 : numel(C)
    r = max(r, size(C{i}, 2));
end

end
