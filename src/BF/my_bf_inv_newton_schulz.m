function [K_inv_BF, history] = ...
    my_bf_inv_newton_schulz(...
    n, K_BF, n_leaf_bf, r_ns, tol_ns, maxit, num_test, ...
    d, type_bf, type_bf_inv, verbose)

arguments (Input)
    n (1, 1) double;
    K_BF;
    n_leaf_bf (1, 1) double;
    r_ns (1, 1) double;
    tol_ns (1, 1) double;
    maxit (1, 1) double;
    num_test (1, 1) double;
    d (1, 1) double;
    type_bf (1, 1) string;
    type_bf_inv (1, 1) string;
    verbose (1, 1) double = 1;
end

arguments (Output)
    K_inv_BF;
    history (1, 1) struct;
end

all_type_bf = ["bf", "mbf", "pbf", "fbf", ...
    "fmbf", "fpbf", "mybf"];
assert(n >= 1 && n == floor(n), ...
    "FIO:BFInverse:InvalidSize", ...
    "n must be a positive integer.");
assert(n_leaf_bf >= 0 && n_leaf_bf == floor(n_leaf_bf), ...
    "FIO:BFInverse:InvalidLeafSize", ...
    "n_leaf_bf must be a nonnegative integer.");
assert(any(type_bf == all_type_bf), ...
    "FIO:BFInverse:UnsupportedInputType", ...
    "Unsupported input BF type: %s.", type_bf);

switch type_bf_inv
    case "bf"
        compatible_type_bf = ["bf", "fbf"];
    case "mbf"
        compatible_type_bf = ["mbf", "fmbf"];
    case "pbf"
        compatible_type_bf = ["pbf", "fpbf"];
    case "mybf"
        compatible_type_bf = "mybf";
    otherwise
        error("FIO:BFInverse:UnsupportedType", ...
            "Unsupported inverse BF type: %s.", type_bf_inv);
end
assert(any(type_bf == compatible_type_bf), ...
    "FIO:BFInverse:IncompatibleTypes", ...
    "type_bf_inv = %s is incompatible with type_bf = %s.", ...
    type_bf_inv, type_bf);

switch type_bf_inv
    case {"bf", "mbf", "pbf"}
        if d == 1
            assert(type_bf_inv == "bf", ...
                "FIO:BFInverse:UnsupportedOutputDimension", ...
                "A 1D SMP inverse must use type_bf_inv = bf.");
            [K_inv_BF, history] = BF_Inv_Newton_Schulz(...
                n, K_BF, r_ns, tol_ns, maxit, num_test, ...
                verbose, type_bf, type_bf_inv);
        elseif d == 2
            [K_inv_BF, history] = BF_Inv_Newton_Schulz2D(...
                n, K_BF, r_ns, tol_ns, maxit, num_test, ...
                verbose, type_bf, type_bf_inv);
        else
            error("FIO:BFInverse:UnsupportedDimension", ...
                "Only 1D and 2D problems are supported.");
        end
    case "mybf"
        if d == 1
            assert(isa(K_BF, "Butterfly"), ...
                "FIO:BFInverse:InvalidButterfly", ...
                "A 1D Butterfly object is required.");
            assert(K_BF.N_ == n, "FIO:BFInverse:SizeMismatch", ...
                "n must agree with the size of the Butterfly object.");
        elseif d == 2
            assert(isa(K_BF, "Butterfly2D"), ...
                "FIO:BFInverse:InvalidButterfly2D", ...
                "A 2D Butterfly object is required.");
            assert(K_BF.n_ == n, "FIO:BFInverse:SizeMismatch", ...
                "n must agree with the size of the Butterfly2D object.");
        else
            error("FIO:BFInverse:UnsupportedDimension", ...
                "Only 1D and 2D problems are supported.");
        end

        [K_inv_BF, history] = K_BF.NewtonSchulzInverse(...
            r_ns, tol_ns, maxit, num_test, verbose);
        history.output_type_bf = "mybf";
end

history.requested_output_type_bf = type_bf_inv;

end
