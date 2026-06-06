classdef Butterfly2D < handle

    properties
        n_ (1, 1) double;
        N_ (1, 1) double;

        num_children_ (1, 1) double;
        ch_list_ (1, :) double;

        tree_ (1, :) cell;
        L_ (1, 1) double;
        h_x_ (1, 1) double;
        h_xi_ (1, 1) double;
        L_x_ (1, 1) double;
        L_xi_ (1, 1) double;

        mid_row_size_ (1, 1) double;
        mid_col_size_ (1, 1) double;

        x_perm_ (:, 1) double;
        x_perm_inv_ (:, 1) double;

        xi_perm_ (:, 1) double;
        xi_perm_inv_ (:, 1) double;
        
        U_ (:, :) cell;
        G_ (:, :) cell;
        M_ (:, :) cell;
        H_ (:, :) cell;
        V_ (:, :) cell;
    end

    methods
        function BF = Butterfly2D(n, a_func, phi_func, ...
                min_points, r, tol)
            BF.n_ = n;
            BF.N_ = n^2;
            BF.ConstructTree(n, min_points);
            BF.Construct(a_func, phi_func, r);
            nnz_before = BF.Nnz();
            BF.OutCompression(tol);
            BF.InCompression(tol);
            nnz_after = BF.Nnz();
            compression_ratio = nnz_before / nnz_after;
            fprintf("  compression_ratio: %.1e\n", compression_ratio);
            BF.SetMidSize();
            fprintf("  n: %d, mid_row_size: %d, mid_col_size: %d\n", ...
                n, BF.mid_row_size_, BF.mid_col_size_);
        end
    end
end