function BBC_ConstructGenerators(A, level, target_rank, tol)
% BBC_ConstructGenerators

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
    target_rank (1, 1) double;
    tol (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 0
        A.BBC_MergeAuxiliaryMatrix();
    end
    target_rank = min(target_rank, size(A.BBC_Y_, 1));
    P = NullBasis(A.BBC_Omega_, target_rank + 5);
    [A.Umat_, A.rank_] = ColBasis(A.BBC_Y_ * P, target_rank, tol);
    Y_OmegaInv = A.BBC_Y_ / A.BBC_Omega_;
    tmp = Y_OmegaInv - A.Umat_ * (A.Umat_' * Y_OmegaInv);
    A.BBC_A_ = tmp + A.Umat_ * (A.Umat_' * tmp');
    if A.leaf_ == 1
        A.Amat_ = A.BBC_A_;
    else
        % Assign R and W.
        offset = 0;
        for i = 1 : A.num_children_
            current_size = A.children_{i}.rank_;
            A.Rmat_{i} = A.Umat_(...
                (offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
        end

        % Assign B.
        A.Bmat_ = cell(A.num_children_, A.num_children_);
        row_offset = 0;
        for i = 1 : A.num_children_
            current_row_size = A.children_{i}.rank_;
            col_offset = 0;
            for j = 1 : i
                current_col_size = A.children_{j}.rank_;
                A.Bmat_{i, j} = A.BBC_A_(...
                    (row_offset + 1) : (row_offset + current_row_size), ...
                    (col_offset + 1) : (col_offset + current_col_size));
                col_offset = col_offset + current_col_size;
            end
            row_offset = row_offset + current_row_size;
        end
    end
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.BBC_ConstructGenerators(level, target_rank, tol);
    end
end

end