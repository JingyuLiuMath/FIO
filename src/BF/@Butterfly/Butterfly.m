classdef Butterfly < handle

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

        constructed_ (1, 1) logical;
        
        U_ (:, :) cell;
        G_ (:, :) cell;
        M_ (:, :) cell;
        H_ (:, :) cell;
        V_ (:, :) cell;
    end

    methods
        function BF = Butterfly(n, n_leaf)
            arguments (Input)
                n (1, 1) double;
                n_leaf (1, 1) double;
            end

            BF.n_ = n;
            BF.N_ = n;
            BF.ConstructTree(n_leaf);
            BF.constructed_ = false;
            BF.U_ = {};
            BF.G_ = {};
            BF.M_ = {};
            BF.H_ = {};
            BF.V_ = {};
        end
    end
end
