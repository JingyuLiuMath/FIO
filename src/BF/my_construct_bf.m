function [BF, Rcomp] = my_construct_bf(k_func, x, xi, r_bf, tol_bf, type_bf)

switch type_bf
    case "bf"
        [BF, Rcomp] = fastBF(k_func, x, xi, r_bf, tol_bf);
    case "mbf"
        [BF, Rcomp] = fastMBF(k_func, x, xi, r_bf, tol_bf);
end

end
