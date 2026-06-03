function run_FIO2D_debug_mbf(result)

n = result.n;
N = result.N;

r_bf = result.r_bf;
tol_bf = result.tol_bf;
type_bf = "mbf";

fprintf("Basic info.\n");
fprintf("  n: %d\n", n);
fprintf("  N: %d\n", N);

fprintf("  r_bf: %d\n", r_bf);
fprintf("  tol_bf: %.1e\n", tol_bf);

fprintf("  t_construct_BF: %.1e\n", result.t_construct_BF);
fprintf("  t_apply_BF: %.1e\n", result.t_apply_BF);
fprintf("  rel_err_BF: %.1e\n", result.rel_err_BF);

M = my_getM_bf(result.K_BF, type_bf);
fprintf("  M in BF: %d\n", M);

nnz_BF = my_nnz_bf(result.K_BF, type_bf);

end