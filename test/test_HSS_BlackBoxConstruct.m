%% Setting.
clear;
close all;

p = 7;
n = 2^p;
N = n^2;
min_points = 64;

tol_hss = 1e-3;

hss_rank_rule = "sqrt";

switch hss_rank_rule
    case "const"
        r_hss = 10;
    case "log"
        r_hss = ceil(log2(N));
    case "sqrt"
        r_hss = ceil(sqrt(N));
end

A = BF_HSS2D(n, n);
[~] = A.BuildTree(min_points);
A.RandInit(hss_rank_rule);
op_A = @(v) A.Apply(v);

B = BF_HSS2D(n, n);
[~] = B.BuildTree(min_points);
B.BlackBoxConstruct(op_A, r_hss, tol_hss);

f = randn(N, 1) + 1i * randn(N, 1);
Af = A.Apply(f);
Bf = B.Apply(f);
rel_err = norm(Af - Bf) / norm(Af);
fprintf("rel_err: %.1e\n", rel_err);

