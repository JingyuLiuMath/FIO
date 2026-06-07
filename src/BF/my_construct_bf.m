function BF = my_construct_bf(n, a_func, phi_func, ...
    x, xi, ...
    n_leaf_bf, r_bf, tol_bf, ...
    type_bf)

exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));
k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

d = size(x, 2);

switch type_bf
    case "bf"
        [BF, ~] = fastBF(k_func, x, xi, r_bf, tol_bf);
    case "mbf"
        [BF, ~] = fastMBF(k_func, x, xi, r_bf, tol_bf);
    case "pbf"
        [BF, ~] = fastBF(k_func, x, xi, r_bf, tol_bf, "polar");
    case "mybf"
        if d == 1
            BF = Butterfly(n, a_func, phi_func, n_leaf_bf, r_bf, tol_bf);
        elseif d == 2
            BF = Butterfly2D(n, a_func, phi_func, n_leaf_bf, r_bf, tol_bf);
        end
end

end
