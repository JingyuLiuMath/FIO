function ULV_Solve_MergeVector(A)
% ULV_Solve_MergeVector

arguments (Input)
    A HSS_Herm;
end

row_size = 0;
for i = 1 : A.num_children_
    row_size = row_size + A.children_{i}.ULV_sk_size_;
end

A.fvec_ = zeros(row_size, A.vec_col_size_);
offset = 0;
for i = 1 : A.num_children_
    current_size = A.children_{i}.ULV_sk_size_;
    A.fvec_((offset + 1) : (offset + current_size), :) ...
        = A.children_{i}.ULV_f_sk_;
    offset = offset + current_size;

    % Clear.
    A.children_{i}.ULV_f_sk_ = [];
end

end