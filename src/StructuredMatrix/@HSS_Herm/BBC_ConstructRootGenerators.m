function BBC_ConstructRootGenerators(A)
% BBC_ConstructRootGenerators

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 1
    A.Amat_ = A.BBC_Y_ / A.BBC_Omega_;
    A.Amat_ = (A.Amat_ + A.Amat_') / 2;
else
    A.BBC_MergeAuxiliaryMatrix();
    A.BBC_A_ = A.BBC_Y_ / A.BBC_Omega_;

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

% Clear.
A.BBC_Omega_ = [];
A.BBC_Y_ = [];
A.BBC_A_ = [];

end