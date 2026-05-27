function res = fun_1D(x, xi)

phi_lin = x * xi.';

sx = (2 + sin(2 * pi * x)) / 8;
phi_nonlin = sx * abs(xi).';

phi = phi_lin + phi_nonlin;

res = complex(cos(2 * pi * phi), sin(2 * pi * phi));

end
