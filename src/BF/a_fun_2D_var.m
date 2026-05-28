function res = a_fun_2D_var(x, xi)

factor = 16;

c1_func = @(x1, x2) (2 + sin(2 * pi * x1) .* sin(2 * pi * x2)) / factor;
c2_func = @(x1, x2) (2 + cos(2 * pi * x1) .* cos(2 * pi * x2)) / factor;

rho_func = @(x1, x2, xi1, xi2) c1_func(x1, x2) * sqrt((xi1.^2)' + (xi2.^2)');

phi_nonlin = rho_func(x(:, 1), x(:, 2), xi(:, 1), xi(:, 2));

phi = phi_nonlin;
[row_ind, col_ind] = find(phi == 0);

res = besselh(0, 2 * pi * phi) .* exp(-2 * pi * 1i * phi);

for it = 1 : length(row_ind)
    j = row_ind(it);
    k = col_ind(it);
    x_j = x(j, :);
    xi_k = xi(k, :);
    my_rho_func = @(xi1, xi2) rho_func(x_j(1), x_j(2), xi1, xi2).';
    my_func = @(xi1, xi2) ...
        besselh(0, 2 * pi * my_rho_func(xi1, xi2)) ...
        .* exp(-2 * pi * my_rho_func(xi1, xi2));
    res(j, k) = integral2(my_func , ...
        xi_k(1) - 0.5, xi_k(1) + 0.5, ...
        xi_k(2) - 0.5, xi_k(2) + 0.5);
end

end