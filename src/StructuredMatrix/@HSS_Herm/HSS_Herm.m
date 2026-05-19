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

        level_size_ (1, 1) double;
        %------------------------------------------------------------------

    end

end