function res = fun_1D(x, xi)

factor = 8;

c1_func = @(x1) (2 + sin(2 * pi * x1)) / factor;

rho_func = @(x1, xi1) c1_func(x1) * abs(xi1).';

phi_lin = x * xi.';

phi_nonlin = rho_func(x, xi);

phi = phi_lin + phi_nonlin;

res = complex(cos(2 * pi * phi), sin(2 * pi * phi));

end
