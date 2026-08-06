function result_list = test_FIO2D_inv_HSS_to_BF(p_list)
% test_FIO2D_inv_HSS_to_BF compresses an HSS approximate inverse into a BF.

arguments (Input)
    p_list (:, 1) double = (4 : 5)';
end

arguments (Output)
    result_list (:, 1) struct;
end

rng(1);

a_func = @(x, xi) ones(size(x, 1), size(xi, 1));
phi_func = @(x, xi) phi_fun_2D(x, xi);

n_leaf_bf = 8;
r_bf = 10;
tol_bf = 1e-8;
type_bf = "fbf";

N_leaf_hss = 64;
tol_hss = 1e-3;
rank_func_tol = @(tol) ceil(3 * log10(1 / tol));

n_leaf_bf_inv = 0;
r_bf_inv = 64;
tol_bf_inv = 1e-3;
type_bf_inv = "bf";

num_test = 4;
result_cell = cell(length(p_list), 1);

for it_p = 1 : length(p_list)
    fprintf("\n\n");
    result_cell{it_p} = RunCase(...
        p_list(it_p), a_func, phi_func, ...
        n_leaf_bf, r_bf, tol_bf, type_bf, ...
        N_leaf_hss, rank_func_tol, tol_hss, ...
        n_leaf_bf_inv, r_bf_inv, tol_bf_inv, type_bf_inv, ...
        num_test);
end

result_list = vertcat(result_cell{:});

fprintf("\nSummary.\n");
fprintf("  n      N   t_HSS(s)  t_ULV(s)  t_BF(s)" ...
    + "  e_BF(op) e_HSS(dir)  e_BF(dir)  nnz(BF)/N^2\n");
for it_p = 1 : length(result_list)
    result = result_list(it_p);
    fprintf("%3d %6d  %8.2e  %8.2e %8.2e" ...
        + "  %8.2e   %8.2e   %8.2e      %8.2e\n", ...
        result.n, result.N, result.t_construct_HSS, ...
        result.t_factor_HSS, result.t_construct_BF_inv, ...
        result.rel_err_BF_inv, result.rel_err_direct_HSS, ...
        result.rel_err_direct_BF, result.ratio_BF_inv);
end

end

function result = RunCase(...
    p, a_func, phi_func, ...
    n_leaf_bf, r_bf, tol_bf, type_bf, ...
    N_leaf_hss, rank_func_tol, tol_hss, ...
    n_leaf_bf_inv, r_bf_inv, tol_bf_inv, type_bf_inv, ...
    num_test)

n = 2^p;
N = n^2;
half_n = n / 2;

fprintf("HSS inverse to BF.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);
fprintf("  r_bf_inv: %d\n", r_bf_inv);
fprintf("  tol_bf_inv: %.1e\n", tol_bf_inv);

x_co = (0 : (n - 1))' / n;
x = TensorProduct2D(x_co, x_co);
xi_co = (-half_n : (half_n - 1))';
xi = TensorProduct2D(xi_co, xi_co);

fprintf("FBF.\n");
t_start = tic;
K_BF = my_construct_bf(...
    n, a_func, phi_func, ...
    x, xi, ...
    n_leaf_bf, r_bf, tol_bf, type_bf);
t_construct_BF = toc(t_start);
fprintf("  t_construct_BF: %.2e\n", t_construct_BF);

op_K = @(v) my_apply_bf(K_BF, v, type_bf);
op_K_adjoint = @(v) my_apply_bf_adj(K_BF, v, type_bf);
op_G = @(v) op_K_adjoint(op_K(v));

fprintf("HSS.\n");
test_G = randn(N, num_test) + 1i * randn(N, num_test);
G_test = op_G(test_G);
c = rank_func_tol(tol_hss);
rel_err_HSS = inf;
num_hss_attempts = 0;

while rel_err_HSS >= 10 * tol_hss
    num_hss_attempts = num_hss_attempts + 1;
    G_HSS = BF_HSS2D(n, n);
    G_HSS.BuildTree(N_leaf_hss);
    rank_func = @(ell) c * n / 2^ell;

    t_start = tic;
    out_HSS = G_HSS.BlackBoxConstruct_Indep_FastBF(...
        op_G, rank_func, tol_hss);
    t_construct_HSS = toc(t_start);

    rel_err_HSS = norm(G_HSS.MyApply(test_G) - G_test, "fro") ...
        / norm(G_test, "fro");
    fprintf("  c: %d\n", c);
    fprintf("  t_construct_HSS: %.2e\n", t_construct_HSS);
    fprintf("  rel_err_HSS: %.2e\n", rel_err_HSS);

    c = 2 * c;
end

hss_rank = G_HSS.Rank();
hss_mem = G_HSS.Storage();
fprintf("  HSS rank: %d\n", hss_rank);
fprintf("  HSS storage / dense: %.2e\n", hss_mem / (16 * N^2));

t_start = tic;
G_HSS.ULV_Factor();
t_factor_HSS = toc(t_start);
fprintf("  t_factor_HSS: %.2e\n", t_factor_HSS);

op_K_inv = @(v) G_HSS.Solve(op_K_adjoint(v));
op_K_inv_adjoint = @(v) op_K(G_HSS.Solve(v));

fprintf("Inverse BF black-box construction.\n");
t_start = tic;
[K_inv_BF, out_BF_inv] = my_construct_bf_blackbox(...
    n, op_K_inv, op_K_inv_adjoint, ...
    x, xi, ...
    n_leaf_bf_inv, r_bf_inv, tol_bf_inv, type_bf_inv);
t_construct_BF_inv = toc(t_start);
fprintf("  t_construct_BF_inv: %.2e\n", t_construct_BF_inv);
fprintf("  operator vectors: %d + %d\n", ...
    out_BF_inv.num_operator_vectors, ...
    out_BF_inv.num_operator_vectors_adj);

test_BF = randn(N, num_test) + 1i * randn(N, num_test);
K_inv_test = op_K_inv(test_BF);
t_start = tic;
K_inv_test_BF = my_apply_bf(K_inv_BF, test_BF, type_bf_inv);
t_apply_BF_inv = toc(t_start);
rel_err_BF_inv = norm(K_inv_test_BF - K_inv_test, "fro") ...
    / norm(K_inv_test, "fro");

K_inv_adjoint_test = op_K_inv_adjoint(test_BF);
K_inv_adjoint_test_BF = my_apply_bf_adj(...
    K_inv_BF, test_BF, type_bf_inv);
rel_err_BF_inv_adjoint = norm(...
    K_inv_adjoint_test_BF - K_inv_adjoint_test, "fro") ...
    / norm(K_inv_adjoint_test, "fro");

f_ex = randn(N, num_test) + 1i * randn(N, num_test);
rhs = op_K(f_ex);
t_start = tic;
f_direct_HSS = op_K_inv(rhs);
t_direct_HSS = toc(t_start);
t_start = tic;
f_direct_BF = my_apply_bf(K_inv_BF, rhs, type_bf_inv);
t_direct_BF = toc(t_start);

rel_err_direct_HSS = norm(f_direct_HSS - f_ex, "fro") ...
    / norm(f_ex, "fro");
rel_err_direct_BF = norm(f_direct_BF - f_ex, "fro") ...
    / norm(f_ex, "fro");
rel_res_direct_BF = norm(op_K(f_direct_BF) - rhs, "fro") ...
    / norm(rhs, "fro");

nnz_BF_inv = my_nnz_bf(K_inv_BF, type_bf_inv);
ratio_BF_inv = nnz_BF_inv / N^2;

fprintf("Checks.\n");
fprintf("  rel_err_BF_inv: %.2e\n", rel_err_BF_inv);
fprintf("  rel_err_BF_inv_adjoint: %.2e\n", ...
    rel_err_BF_inv_adjoint);
fprintf("  rel_err_direct_HSS: %.2e\n", rel_err_direct_HSS);
fprintf("  rel_err_direct_BF: %.2e\n", rel_err_direct_BF);
fprintf("  rel_res_direct_BF: %.2e\n", rel_res_direct_BF);
fprintf("  nnz(BF_inv) / N^2: %.2e\n", ratio_BF_inv);
fprintf("  t_direct_HSS: %.2e\n", t_direct_HSS);
fprintf("  t_direct_BF: %.2e\n", t_direct_BF);

result = struct();
result.p = p;
result.n = n;
result.N = N;
result.r_bf_inv = r_bf_inv;
result.tol_bf_inv = tol_bf_inv;
result.t_construct_BF = t_construct_BF;
result.t_construct_HSS = t_construct_HSS;
result.t_factor_HSS = t_factor_HSS;
result.t_construct_BF_inv = t_construct_BF_inv;
result.t_apply_BF_inv = t_apply_BF_inv;
result.t_direct_HSS = t_direct_HSS;
result.t_direct_BF = t_direct_BF;
result.rel_err_HSS = rel_err_HSS;
result.rel_err_BF_inv = rel_err_BF_inv;
result.rel_err_BF_inv_adjoint = rel_err_BF_inv_adjoint;
result.rel_err_direct_HSS = rel_err_direct_HSS;
result.rel_err_direct_BF = rel_err_direct_BF;
result.rel_res_direct_BF = rel_res_direct_BF;
result.hss_rank = hss_rank;
result.hss_mem = hss_mem;
result.nnz_BF_inv = nnz_BF_inv;
result.ratio_BF_inv = ratio_BF_inv;
result.num_hss_attempts = num_hss_attempts;
result.out_HSS = out_HSS;
result.out_BF_inv = out_BF_inv;

end
