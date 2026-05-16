function U_level_list = BBC_Indep_ConstructGenerators(...
    A, level, ...
    Omega, Y, ...
    target_rank, ...
    U_level_list, ...
    tol)
% BBC_ConstructGenerators

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
    Omega (:, :) double;
    Y (:, :) double;
    target_rank (1, 1) double;
    U_level_list (1, :) cell;
    tol (1, 1) double;
end

arguments (Output)
    U_level_list (1, :) cell;
end

if A.level_ == level
    target_rank = min(target_rank, size(Y, 1));
    P = NullBasis(Omega, target_rank + 5);
    [U, A.rank_] = ColBasis(Y * P, target_rank, tol);
    U_level_list{end + 1} = U;
    Y_OmegaInv = Y / Omega;
    tmp = Y_OmegaInv - U * (U' * Y_OmegaInv);
    Acheck = tmp + U * (U' * tmp');

    if A.leaf_ == 1
        A.Amat_ = Acheck;
        A.Umat_ = U;
    else
        % Assign R and W.
        offset = 0;
        for i = 1 : A.num_children_
            current_size = A.children_{i}.rank_;
            A.Rmat_{i} = U(...
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
        U_level_list = A.children_{i}.BBC_Indep_ConstructGenerators(...
            level, ...
            Omega((offset + 1) : (offset + current_size), :), ...
            Y((offset + 1) : (offset + current_size), :), ...
            target_rank, ...
            U_level_list, ...
            tol);
        offset = offset + current_size;
    end
end

end