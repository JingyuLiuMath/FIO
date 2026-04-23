classdef (Abstract) HSS_Herm < handle
    % HSS

    properties
        % *****************************************************************
        % PROPERTY: Tree.
        level_ (1, 1) double = 0;
        leaf_ (1, 1) double = 0;
        max_level_ (1, 1) double = 0;
        children_ (1, :) cell;
        num_children_ (1, 1) double = 0;
        % -----------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Matrix.
        global_size_ (1, 1) double = 0;
        size_ (1, 1) double = 0;
        rank_ (1, 1) double = 0;  % size(U, 2), rank of the col space.
        offset_ (1, 1) double = 0;
        % -----------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Generators.
        % Basis matrices.
        Umat_ (:, :) double;

        % Transfer matrices.
        Rmat_ (1, :) cell;

        % Interaction matrices.
        Bmat_ (:, :) cell;

        % Full matrices.
        Amat_ (:, :) double;
        %------------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Vectors.
        uvec_ (:, :) double;
        uhvec_ (:, :) double;
        fhvec_ (:, :) double;
        fvec_ (:, :) double;
        vec_col_size_ (:, :) double;
        %------------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: ULV.
        % Matrices.
        ULV_Q_ (:, :) double;  % Col zeroing-out matrix.
        ULV_U_sk_ (:, :) double;
        ULV_re_size_ (1, 1) double;
        ULV_sk_size_ (1, 1) double;
        ULV_A_re_re_ (:, :) double;
        ULV_A_sk_re_ (:, :) double;
        ULV_A_sk_sk_ (:, :) double;
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Vectors.
        ULV_f_sk_ (:, :) double;
        ULV_u_sk_ (:, :) double;
        ULV_u_re_ (:, :) double;
        %------------------------------------------------------------------

        % *****************************************************************
        % PROPERTY: Black-Box Construction.
        BBC_Omega_ (:, :) double;
        BBC_Y_ (:, :) double;
        BBC_A_ (:, :) double;
        %------------------------------------------------------------------

    end

    methods
        % ****************************************************************
        % METHOD: Utilization
        U = Matrix_U(A);
        mem = Storage(A);
        r = Rank(A);
        FillVector_Col(A, u);
        f = FetchVector_Row(A);
        FillVector_Row(A, f);
        u = FetchVector_Col(A);
        max_leaf_size = MaxLeafSize(A);
        %-----------------------------------------------------------------

        % ****************************************************************
        % METHOD: Apply.
        f = Apply(A, u)
        Apply_Upward(A, level);
        Apply_Root(A);
        Apply_Downward(A, level);
        %-----------------------------------------------------------------

        % ****************************************************************
        % METHOD: ULV factorization and the least squares solver.
        % Factorization.
        ULV_Factor(A);
        ULV_Eliminate(A, level);
        ULV_Merge(A);
        ULV_RootFactor(A);
        % Solution.
        u = ULV_Solve(A, f);
        ULV_Solve_Upward(A, level);
        ULV_Solve_Merge(A);
        ULV_Solve_Root(A);
        ULV_Solve_Split(A);
        ULV_Solve_Downward(A, level);
        %-----------------------------------------------------------------

        % ****************************************************************
        % METHOD: Block-Box Construction.
        BlackBoxConstruct(A, op_A, op_Astar, ...
            row_leaf_size, col_leaf_size, target_rank);
        BBC_FillAuxiliaryMatrix(A, Omega, Y, Psi, Z);
        BBC_ConstructGenerators(A, level, target_rank, tol);
        BBC_MergeAuxiliaryMatrix(A);
        BBC_ConstructRootGenerators(A);
        BBC_EliminateRootBMat(A);
        BBC_EliminateBMat(A, level);
        %-----------------------------------------------------------------

    end

end