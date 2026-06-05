classdef Butterfly < handle

    properties
        N_ (1, 1) double;

        tree_ (1, :) cell;
        L_ (1, 1) double;
        h_x_ (1, 1) double;
        h_xi_ (1, 1) double;
        L_x_ (1, 1) double;
        L_xi_ (1, 1) double;
        
        U_ (:, :) cell;
        G_ (:, :) cell;
        M_ (:, :) cell;
        H_ (:, :) cell;
        V_ (:, :) cell;
    end

    methods
        function BF = Butterfly(N,a_func, phi_func, r_bf, tol_bf)
            BF.N_ = N;
            BF.ConstructTree(N);
            BF.Construct(a_func, phi_func, r_bf);
            nnz_before = BF.Nnz();
            BF.BFOutCompression(tol_bf);
            nnz_after = BF.Nnz();
            compression_ratio = nnz_before / nnz_after;
            fprintf("  compression ratio: %.1e\n", compression_ratio);
        end
    end
end