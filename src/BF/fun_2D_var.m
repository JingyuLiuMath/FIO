function res = fun_2D_var(x, xi)

factor = 16;

phi_lin = x(:, 1) * xi(:, 1)' + x(:, 2) * xi(:, 2)';

sx = (2 + sin(2 * pi * x(:, 1)) .* sin(2 * pi * x(:, 2))) / factor;
rk = sqrt(xi(:, 1).^2 + xi(:, 2).^2);
phi_nonlin = sx * rk.';

phi = phi_lin + phi_nonlin;

res = complex(cos(2 * pi * phi), sin(2 * pi * phi));

end
