%% Setting.
clear;
close all;

p = 6;
n = 2^p;
N = n^2;
min_points = 256;

tol_hss = 1e-14;

hss_rank_rule_list = ["const"; "log"; "sqrt"];
hss_rank_rule = "sqrt";

A = BF_HSS2D(n, n);
A.BuildTree(min_points);
A.RandInit(hss_rank_rule);
fprintf("hss rank: %d\n", A.Rank());
fprintf("max level: %d\n", A.max_level_);
op_A = @(v) A.Apply(v);

c_ex = complex(randn(N, 1), randn(N, 1));
f_ex = op_A(c_ex);

%%
fprintf("\n\n");
fprintf("old\n");
B = BF_HSS2D(n, n);
B.BuildTree(min_points);
rel_err = inf;
r_hss = 10;
tic;
while rel_err > tol_hss * 10
    fprintf("r_hss: %d\n", r_hss);
    B.BlackBoxConstruct(op_A, r_hss, tol_hss);
    f = B.Apply(c_ex);
    rel_err = norm(f - f_ex) / norm(f_ex);
    fprintf("rel_err: %.1e\n", rel_err);
    r_hss = r_hss * 2;
end
t_old = toc;
fprintf("t_old: %.1e\n", t_old);

%%
fprintf("\n\n");
fprintf("new\n");
tic;
rank_func = @(ell) n / 2^ell;
C = BF_HSS2D(n, n);
C.BuildTree(min_points);
C.BlackBoxConstruct_Indep(op_A, rank_func, tol_hss);
f = C.Apply(c_ex);
rel_err = norm(f - f_ex) / norm(f_ex);
fprintf("rel_err: %.1e\n", rel_err);
t_new = toc;
fprintf("t_new: %.1e\n", t_new);