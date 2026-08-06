function [K_inv_BF, history] = BF_Inv_Newton_Schulz2D(...
    n, K_BF, r_ns, tol_ns, maxit, num_test, verbose, ...
    type_bf, type_bf_inv)
% BF_Inv_Newton_Schulz2D constructs a 2D SMP approximate inverse.

arguments (Input)
    n (1, 1) double;
    K_BF;
    r_ns (1, 1) double;
    tol_ns (1, 1) double;
    maxit (1, 1) double;
    num_test (1, 1) double;
    verbose (1, 1) double;
    type_bf (1, 1) string;
    type_bf_inv (1, 1) string;
end

arguments (Output)
    K_inv_BF;
    history (1, 1) struct;
end

r_ns = ceil(r_ns);
maxit = ceil(maxit);
num_test = ceil(num_test);
assert(r_ns >= 1, "FIO:BFInverse2D:InvalidRank", ...
    "r_ns must be positive.");
assert(tol_ns >= 0, "FIO:BFInverse2D:InvalidNSTolerance", ...
    "tol_ns must be nonnegative.");
assert(maxit >= 0, "FIO:BFInverse2D:InvalidMaxit", ...
    "maxit must be nonnegative.");
assert(num_test >= 1, "FIO:BFInverse2D:InvalidNumTest", ...
    "num_test must be positive.");
assert(n >= 1 && n == floor(n), "FIO:BFInverse2D:InvalidSize", ...
    "n must be a positive integer.");
assert(any(type_bf == ["bf", "mbf", "pbf", "fbf", ...
    "fmbf", "fpbf", "mybf"]), ...
    "FIO:BFInverse2D:UnsupportedType", ...
    "Unsupported 2D input BF type: %s.", type_bf);
assert(any(type_bf_inv == ["bf", "mbf", "pbf"]), ...
    "FIO:BFInverse2D:UnsupportedOutputType", ...
    "type_bf_inv must be bf, mbf, or pbf for a 2D SMP inverse.");

if type_bf == "mybf"
    assert(isa(K_BF, "Butterfly2D"), ...
        "FIO:BFInverse2D:InvalidButterfly2D", ...
        "A 2D Butterfly object is required for type_bf = mybf.");
    assert(K_BF.n_ == n, "FIO:BFInverse2D:SizeMismatch", ...
        "n must agree with the size of the Butterfly2D object.");
end
N = n^2;
half_n = n / 2;
x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);
xi_co = (-half_n : (half_n - 1))';
xi = TensorProduct2D(xi_co, xi_co);

num_power_iter = 10;
alpha_safety = 0.95;
t_total_start = tic;
t_power_start = tic;

q = randn(N, 1);
q = q / norm(q);
for it = 1 : num_power_iter
    z = ApplyKAdjoint(...
        K_BF, ApplyK(K_BF, q, type_bf), type_bf);
    norm_z = norm(z);
    assert(norm_z > 0, "FIO:BFInverse2D:SingularOperator", ...
        "The operator appears to be zero.");
    q = z / norm_z;
end
sigma_max_est = norm(ApplyK(K_BF, q, type_bf));
alpha = alpha_safety / sigma_max_est^2;
t_power = toc(t_power_start);

if verbose == 1
    fprintf("  Newton-Schulz 2D SMP inverse.\n");
    fprintf("    estimated norm: %.3e\n", sigma_max_est);
    fprintf("    initial scaling: %.3e\n", alpha);
end

t_construct_start = tic;
K_inv_BF = AdjointView(K_BF, type_bf, alpha);
type_X = "adjoint_view";
t_construct = toc(t_construct_start);

test_mat = randn(N, num_test);
test_norm = norm(test_mat, "fro");
K_test_mat = ApplyK(K_BF, test_mat, type_bf);

history.iteration = (0 : maxit)';
history.residual = nan(maxit + 1, 1);
history.nnz = nan(maxit + 1, 1);
history.construct_time = nan(maxit + 1, 1);
history.compression_tolerance = nan(maxit + 1, 1);
history.sigma_max_est = sigma_max_est;
history.alpha = alpha;
history.power_iteration_time = t_power;
history.input_type_bf = type_bf;
history.requested_output_type_bf = type_bf_inv;
history.initialization_method = "scaled adjoint view";
history.initialization_time = t_construct;
history.residual_time = 0;
history.blackbox_stats = cell(maxit + 1, 1);

t_residual_start = tic;
residual = Residual(...
    K_inv_BF, test_mat, K_test_mat, test_norm, type_X);
history.residual_time = history.residual_time + toc(t_residual_start);
history.residual(1) = residual;
history.nnz(1) = my_nnz_bf(K_inv_BF, type_X);
history.construct_time(1) = t_construct;
history.compression_tolerance(1) = nan;
history.blackbox_stats{1} = [];

if verbose == 1
    fprintf("    iter %2d: residual %.3e, " + ...
        "tol %.3e, ratio %.3e\n", ...
        0, residual, history.compression_tolerance(1), ...
        history.nnz(1) / N^2);
end

last_iter = 0;
for iter = 1 : maxit
    if residual <= tol_ns
        break;
    end

    X_old = K_inv_BF;
    type_X_old = type_X;
    op_X_next = @(f) NewtonSchulzUpdate(...
        K_BF, X_old, f, type_bf, type_X_old);
    op_X_next_adjoint = @(f) ...
        NewtonSchulzUpdateAdjoint(...
        K_BF, X_old, f, type_bf, type_X_old);

    residual_previous = residual;
    tol_compression = max(tol_ns, min(...
        2 * tol_ns, 2 * residual_previous^2));

    t_construct_start = tic;
    [K_inv_BF, out_blackbox] = ConstructInverseBF(...
        n, op_X_next, op_X_next_adjoint, ...
        x, xi, r_ns, tol_compression, type_bf_inv);
    type_X = type_bf_inv;
    t_construct = toc(t_construct_start);

    t_residual_start = tic;
    residual = Residual(...
        K_inv_BF, test_mat, K_test_mat, test_norm, type_X);
    history.residual_time = history.residual_time ...
        + toc(t_residual_start);
    history.residual(iter + 1) = residual;
    history.nnz(iter + 1) = ...
        my_nnz_bf(K_inv_BF, type_X);
    history.construct_time(iter + 1) = t_construct;
    history.compression_tolerance(iter + 1) = tol_compression;
    history.blackbox_stats{iter + 1} = out_blackbox;
    last_iter = iter;

    if verbose == 1
        fprintf("    iter %2d: residual %.3e, " + ...
            "tol %.3e, ratio %.3e\n", ...
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
history.blackbox_stats = history.blackbox_stats(1 : last_ind);
history.num_iterations = last_iter;
history.num_blackbox_constructions = last_iter;
history.converged = residual <= tol_ns;
history.output_type_bf = type_X;
history.total_time = toc(t_total_start);
history.nonconstruction_time = max(...
    0, history.total_time - sum(history.construct_time));

if verbose == 1
    fprintf("    converged: %d\n", history.converged);
end

if type_X == "adjoint_view"
    if history.converged
        warning("FIO:BFInverse2D:InitialGuessRetained", ...
            "The scaled adjoint already satisfies tol_ns. " + ...
            "It is returned as an adjoint view without a " + ...
            "%s black-box construction.", type_bf_inv);
    elseif maxit == 0
        warning("FIO:BFInverse2D:ZeroIterationStructureRetained", ...
            "maxit is zero. The scaled adjoint is returned as an " + ...
            "adjoint view without a %s black-box construction.", ...
            type_bf_inv);
    end
end

end

function [Factor, out] = ConstructInverseBF(...
    n, op_K, op_K_adjoint, x, xi, r_bf, tol_bf, type_bf_inv)

[Factor, out] = my_construct_bf_blackbox(...
    n, op_K, op_K_adjoint, ...
    x, xi, 0, r_bf, tol_bf, type_bf_inv);

end

function X_BF = AdjointView(K_BF, type_bf, alpha)

X_BF = struct();
X_BF.factor = K_BF;
X_BF.type_bf = type_bf;
X_BF.alpha = alpha;

end

function y = ApplyK(K_BF, y, type_bf)

y = my_apply_bf(K_BF, y, type_bf);

end

function y = ApplyKAdjoint(K_BF, y, type_bf)

y = my_apply_bf_adj(K_BF, y, type_bf);

end

function y = NewtonSchulzUpdate(...
    K_BF, X_BF, y, type_bf, type_X)

u = ApplyX(X_BF, y, type_X);
y = 2 * u - ApplyX(...
    X_BF, ApplyK(K_BF, u, type_bf), type_X);

end

function y = NewtonSchulzUpdateAdjoint(...
    K_BF, X_BF, y, type_bf, type_X)

u = ApplyXAdjoint(X_BF, y, type_X);
y = 2 * u - ApplyXAdjoint(...
    X_BF, ApplyKAdjoint(K_BF, u, type_bf), type_X);

end

function residual = Residual(...
    X_BF, test_mat, K_test_mat, test_norm, type_X)

residual = norm(...
    ApplyX(X_BF, K_test_mat, type_X) ...
    - test_mat, "fro") ...
    / test_norm;

end

function y = ApplyX(X_BF, y, type_X)

y = my_apply_bf(X_BF, y, type_X);

end

function y = ApplyXAdjoint(X_BF, y, type_X)

y = my_apply_bf_adj(X_BF, y, type_X);

end
