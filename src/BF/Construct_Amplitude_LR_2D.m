function [U, V] = Construct_Amplitude_LR_2D(a_func, x, xi, r_est, tol)

N = size(x, 1);

k_est = 3 * r_est;

col_index = randperm(N, k_est);
A_J_C = a_func(x, xi(col_index, :));
[U, ~ ,~] = MySVDSketch(A_J_C, tol);

row_index = randperm(N, k_est);
A_R_K = a_func(x(row_index, :), xi);
U_R = U(row_index, :);

V = (U_R \ A_R_K)';

end