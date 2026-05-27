function res = a_fun_2D_var(x, xi)

factor = 16;

sx = (2 + sin(2 * pi * x(:, 1)) .* sin(2 * pi * x(:, 2))) / factor;
rk = sqrt(xi(:, 1).^2 + xi(:, 2).^2);
phi_nonlin = sx * rk.';

phi = phi_nonlin;
[row_ind, col_ind] = find(phi == 0);

res = besselh(0, 2 * pi * phi) .* exp(-2 * pi * 1i * phi);

for it = 1 : length(row_ind)
    j = row_ind(it);
    k = col_ind(it);
    x_j = x(j, :);
    xi_k = xi(k, :);
    sx_j = (3 + sin(2 * pi * x_j(1)) .* sin(2 * pi * x_j(2))) / factor;
    rho_func = @(xi1, xi2) sx_j * sqrt(xi1.^2 + xi2.^2);
    a_func_x_j = @(xi1, xi2) ...
        besselh(0, 2 * pi * rho_func(xi1, xi2)) .* exp(-2 * pi * rho_func(xi1, xi2));
    res(j, k) = integral2(a_func_x_j , ...
        xi_k(1) - 0.5, xi_k(1) + 0.5, ...
        xi_k(2) - 0.5, xi_k(2) + 0.5);
end

end