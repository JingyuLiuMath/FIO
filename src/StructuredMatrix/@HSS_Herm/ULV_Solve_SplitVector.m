function ULV_Solve_SplitVector(A)
% ULV_Solve_SplitVector

% Jingyu Liu, November 23, 2024.

arguments (Input)
    A HSS_Herm;
end

offset = 0;
for i = 1 : A.num_children_
    current_size = A.children_{i}.ULV_sk_size_;
    A.children_{i}.ULV_u_sk_ ...
        = A.uvec_((offset + 1) : (offset + current_size), :);
    offset = offset + current_size;
end

% Clear.
A.uvec_ = [];

end