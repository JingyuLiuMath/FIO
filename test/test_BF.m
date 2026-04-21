%% Setting.
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

%% BF
K_BF = bf_explicit(BF_input.fun, ...
    BF_input.xx, BF_input.xbox, ...
    BF_input.kk, BF_input.kbox, ...
    BF_input.mR, BF_input.tol, 0);

%% Apply.
f = randn(N,1) + 1i*randn(N,1);
Kf_ex = K * f;
Kf = apply_bf(K_BF, f);
e = Kf_ex - Kf;
rel_err = norm(e) / norm(Kf_ex);
fprintf("rel err: %.1e\n", rel_err);

f = randn(N,1) + 1i*randn(N,1);
Khf_ex = K' * f;
Khf = apply_bf_adj(K_BF, f);
e = Khf_ex - Khf;
rel_err = norm(e) / norm(Khf_ex);
fprintf("rel err: %.1e\n", rel_err);
