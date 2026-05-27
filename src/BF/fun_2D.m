function res = fun_2D(x, xi)

factor = 16;

c1_func = @(x1, x2) (2 + sin(2 * pi * x1) .* sin(2 * pi * x2)) / factor;
c2_func = @(x1, x2) (2 + cos(2 * pi * x1) .* cos(2 * pi * x2)) / factor;

rho_func = @(x1, x2, xi1, xi2) sqrt(...
    c1_func(x1, x2).^2 * (xi1.^2)' ...
    + c2_func(x1, x2).^2 * (xi2.^2)');

phi_lin = x(:, 1) * xi(:, 1)' + x(:, 2) * xi(:, 2)';

phi_nonlin = rho_func(x(:, 1), x(:, 2), xi(:, 1), xi(:, 2));

phi = phi_lin + phi_nonlin;

res = complex(cos(2 * pi * phi), sin(2 * pi * phi));

end
