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
                n_leaf, r, tol)
            arguments (Input)
                n (1, 1) double;
                a_func function_handle;
                phi_func function_handle;
                n_leaf (1, 1) double;
                r (1, 1) double;
                tol (1, 1) double;
            end

            BF.n_ = n;
            BF.N_ = n^2;
            BF.ConstructTree(n, n_leaf);
            BF.Construct(a_func, phi_func, r);
            BF.OutCompression(tol);
            BF.InCompression(tol);
        end
    end
end