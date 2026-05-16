function BBC_Indep_ConstructGenerators_New(...
    A, level, ...
    Omega, Y, ...
    target_rank, ...
    tol)
% BBC_ConstructGenerators

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
    Omega (:, :) double;
    Y (:, :) double;
    target_rank (1, 1) double;
    tol (1, 1) double;
end

if A.level_ == level
    target_rank = min(target_rank, size(Y, 1));
    P = NullBasis(Omega, target_rank + 5);
    [A.Umat_, A.rank_] = ColBasis(Y * P, target_rank, tol);
    Y_OmegaInv = Y / Omega;
    tmp = Y_OmegaInv - A.Umat_ * (A.Umat_' * Y_OmegaInv);
    Acheck = tmp + A.Umat_ * (A.Umat_' * tmp');
    A.level_size_ = A.rank_;

    if A.leaf_ == 1
        A.Amat_ = Acheck;
    else
        % Assign B.
        A.Bmat_ = cell(A.num_children_, A.num_children_);
        row_offset = 0;
        for i = 1 : A.num_children_
            current_row_size = A.children_{i}.rank_;
            col_offset = 0;
            for j = 1 : i
                current_col_size = A.children_{j}.rank_;
                A.Bmat_{i, j} = Acheck(...
                    (row_offset + 1) : (row_offset + current_row_size), ...
                    (col_offset + 1) : (col_offset + current_col_size));
                col_offset = col_offset + current_col_size;
            end
            row_offset = row_offset + current_row_size;
        end
    end
elseif A.leaf_ == 0
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.level_size_;
        A.children_{i}.BBC_Indep_ConstructGenerators_New(...
            level, ...
            Omega((offset + 1) : (offset + current_size), :), ...
            Y((offset + 1) : (offset + current_size), :), ...
            target_rank, ...
            tol);
        offset = offset + current_size;
    end
end

end