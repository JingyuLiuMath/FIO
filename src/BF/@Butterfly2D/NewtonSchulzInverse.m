function [K_inv_BF, history] = NewtonSchulzInverse(BF, ...
        r_ns, tol_ns, maxit, num_test, verbose)
% NewtonSchulzInverse constructs a BHP approximate inverse.

arguments (Input)
    BF Butterfly2D;
    r_ns (1, 1) double;
    tol_ns (1, 1) double;
    maxit (1, 1) double = 10;
    num_test (1, 1) double = 8;
    verbose (1, 1) double = 1;
end

arguments (Output)
    K_inv_BF Butterfly2D;
    history (1, 1) struct;
end

r_ns = ceil(r_ns);
maxit = ceil(maxit);
num_test = ceil(num_test);
assert(r_ns >= 1, "Butterfly2D:NewtonSchulzInverse:InvalidRank", ...
    "r_ns must be positive.");
assert(tol_ns >= 0, "Butterfly2D:NewtonSchulzInverse:InvalidNSTolerance", ...
    "tol_ns must be nonnegative.");
assert(maxit >= 0, "Butterfly2D:NewtonSchulzInverse:InvalidMaxit", ...
    "maxit must be nonnegative.");
assert(num_test >= 1, "Butterfly2D:NewtonSchulzInverse:InvalidNumTest", ...
    "num_test must be positive.");

N = BF.N_;
num_power_iter = 10;
alpha_safety = 0.95;
t_total_start = tic;
t_power_start = tic;

q = randn(N, 1);
q = q / norm(q);
for it = 1 : num_power_iter
    z = BF.ApplyAdjoint(BF.Apply(q));
    norm_z = norm(z);
    assert(norm_z > 0, "Butterfly2D:NewtonSchulzInverse:SingularOperator", ...
        "The operator appears to be zero.");
    q = z / norm_z;
end
sigma_max_est = norm(BF.Apply(q));
alpha = alpha_safety / sigma_max_est^2;
t_power = toc(t_power_start);

if verbose == 1
    fprintf("  Newton-Schulz BHP inverse.\n");
    fprintf("    estimated norm: %.3e\n", sigma_max_est);
    fprintf("    initial scaling: %.3e\n", alpha);
end

t_construct_start = tic;
K_inv_BF = BF.Adjoint(alpha);
t_construct = toc(t_construct_start);

test_mat = randn(N, num_test);
test_norm = norm(test_mat, "fro");
BF_test_mat = BF.Apply(test_mat);

history.iteration = (0 : maxit)';
history.residual = nan(maxit + 1, 1);
history.nnz = nan(maxit + 1, 1);
history.construct_time = nan(maxit + 1, 1);
history.compression_tolerance = nan(maxit + 1, 1);
history.construct_tolerance = nan(maxit + 1, 1);
history.blackbox_stats = cell(maxit + 1, 1);
history.sigma_max_est = sigma_max_est;
history.alpha = alpha;
history.power_iteration_time = t_power;
history.initialization_method = "scaled adjoint";
history.initialization_time = t_construct;
history.residual_time = 0;

t_residual_start = tic;
residual = Residual(K_inv_BF, test_mat, BF_test_mat, test_norm);
history.residual_time = history.residual_time + toc(t_residual_start);
history.residual(1) = residual;
history.nnz(1) = K_inv_BF.Nnz();
history.construct_time(1) = t_construct;
history.blackbox_stats{1} = [];

if verbose == 1
    fprintf("    iter %2d: residual %.3e, tol %.3e, ratio %.3e\n", ...
        0, residual, history.compression_tolerance(1), ...
        history.nnz(1) / N^2);
end

last_iter = 0;
for iter = 1 : maxit
    if residual <= tol_ns
        break;
    end

    X_old = K_inv_BF;
    op_X_next = @(f) NewtonSchulzUpdate(BF, X_old, f);
    op_X_next_adjoint = @(f) ...
        NewtonSchulzUpdateAdjoint(BF, X_old, f);

    tol_compression = max(tol_ns, min( ...
        2 * tol_ns, 2 * residual^2));
    tol_construct = tol_compression;

    X_new = EmptyButterfly(BF);
    t_construct_start = tic;
    out_blackbox = X_new.BlackBoxConstruct( ...
        op_X_next, op_X_next_adjoint, r_ns, tol_construct, 0);
    t_construct = toc(t_construct_start);
    K_inv_BF = X_new;

    t_residual_start = tic;
    residual = Residual(K_inv_BF, test_mat, BF_test_mat, test_norm);
    history.residual_time = history.residual_time + toc(t_residual_start);
    history.residual(iter + 1) = residual;
    history.nnz(iter + 1) = K_inv_BF.Nnz();
    history.construct_time(iter + 1) = t_construct;
    history.compression_tolerance(iter + 1) = tol_compression;
    history.construct_tolerance(iter + 1) = tol_construct;
    history.blackbox_stats{iter + 1} = out_blackbox;
    last_iter = iter;

    if verbose == 1
        fprintf("    iter %2d: residual %.3e, tol %.3e, ratio %.3e\n", ...
            iter, residual, tol_compression, ...
            history.nnz(iter + 1) / N^2);
    end

    if ~isfinite(residual)
        break;
    end
end

last_ind = last_iter + 1;
history.iteration = history.iteration(1 : last_ind);
history.residual = history.residual(1 : last_ind);
history.nnz = history.nnz(1 : last_ind);
history.construct_time = history.construct_time(1 : last_ind);
history.compression_tolerance = ...
    history.compression_tolerance(1 : last_ind);
history.construct_tolerance = ...
    history.construct_tolerance(1 : last_ind);
history.blackbox_stats = history.blackbox_stats(1 : last_ind);
history.blackbox_total_time = BlackBoxField( ...
    history.blackbox_stats, "total_time");
history.num_iterations = last_iter;
history.num_blackbox_constructions = last_iter;
history.converged = residual <= tol_ns;
history.total_time = toc(t_total_start);
history.non_blackbox_time = max(0, history.total_time ...
    - sum(history.blackbox_total_time));

if verbose == 1
    fprintf("    converged: %d\n", history.converged);
end

end

function value = BlackBoxField(stats_list, field_name)

value = zeros(length(stats_list), 1);
for i = 1 : length(stats_list)
    if ~isempty(stats_list{i})
        value(i) = stats_list{i}.(field_name);
    end
end

end

function y = NewtonSchulzUpdate(BF, X_BF, f)

u = X_BF.Apply(f);
y = 2 * u - X_BF.Apply(BF.Apply(u));

end

function y = NewtonSchulzUpdateAdjoint(BF, X_BF, f)

u = X_BF.ApplyAdjoint(f);
y = 2 * u - X_BF.ApplyAdjoint(BF.ApplyAdjoint(u));

end

function X_BF = EmptyButterfly(BF)

X_BF = Butterfly2D(BF.n_, BF.n_leaf_, 0);

end

function residual = Residual( ...
        X_BF, test_mat, BF_test_mat, test_norm)

residual = norm(X_BF.Apply(BF_test_mat) - test_mat, "fro") ...
    / test_norm;

end
