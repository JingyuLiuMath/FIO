clear;
close all;

c_func = @(x) (2 + sin(2 * pi * x)) / 8;

phi_func = @(x, xi) x * xi.' + c_func(x) * abs(xi).';

N = 2^10;
half_N = N / 2;

func_name = 'fun0';

BF_input = struct();
BF_input.xbox = [1,N+1];
BF_input.xx = (1:N)';

BF_input.kbox = [1,N+1];
BF_input.kk = (1:N)';
switch func_name
    case 'fun0'
        BF_input.fun = @(x,k)fun0(N,x,k);
    case 'funF' % Fourier transform
        BF_input.fun = @(x,k)funF(N,x,k);
    case 'funH' % Hankel function
        BF_input.fun = @(x,k)funH(N,x,k);
end
BF_input.mR = 4;
BF_input.tol = 1e-8;

x = (0 : (N - 1))' / N;
xi = (-half_N : (half_N - 1))';

min_points = 256;
target_rank = 100;
tol = 1e-6;

%% K
K = exp(2 * pi * 1i * phi_func(x, xi));
r_K = rank(K);
fprintf("rank(K): %d\n", r_K);
fprintf("cond(K): %.1e\n", cond(K));
if r_K ~= N
    keyboard;
end

G = K' * K;

%% BF
K_BF = bf_explicit(BF_input.fun, ...
    BF_input.xx, BF_input.xbox, ...
    BF_input.kk, BF_input.kbox, ...
    BF_input.mR, BF_input.tol, 0);

%% BF Apply.
f = randn(N,1) + 1i*randn(N,1);
Kf_ex = K * f;
Kf = apply_bf(K_BF, f);
e = Kf_ex - Kf;
rel_err = norm(e) / norm(Kf_ex);
fprintf("rel err: %.1e\n", rel_err);

%% HSS.
% op_G = @(v) G * v;
op_G = @(v) apply_bf_adj(K_BF, apply_bf(K_BF, v));

G_HSS = BF_HSS(N);
G_HSS.BuildTree(min_points);
leaf_size = G_HSS.MaxLeafSize();
G_HSS.BlackBoxConstruct(op_G, leaf_size, target_rank, tol);
r_hss = G_HSS.Rank();
fprintf("r_hss: %d\n", r_hss);
mem = G_HSS.Storage();
ratio = mem / N^2;
fprintf("ratio: %.1e\n", ratio);

%% ULV.
G_HSS.ULV_Factor();

%%
c_ex = randn(N, 1);

% f = K * c_ex;
% c = G_HSS.ULV_Solve(K' * f);
% res = f - K * c;

f = apply_bf(K_BF, c_ex);
c = G_HSS.ULV_Solve(apply_bf_adj(K_BF, f));
res = f - apply_bf(K_BF, c);

rel_res = norm(res) / norm(f);
fprintf("rel res: %.1e\n", rel_res);
err = c_ex - c;
rel_err = norm(err) / norm(c_ex);
fprintf("rel err: %.1e\n", rel_err);
