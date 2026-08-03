function FillVector_Col(A, u)
% FillVector_Col

arguments (Input)
    A HSS_Herm;
    u (:, :) double;
end

A.vec_col_size_ = size(u, 2);

if A.leaf_ == 1
    A.uvec_ = u;
else
    % Recursion.
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.size_;
        A.children_{i}.FillVector_Col(...
            u((offset + 1) : (offset + current_size), :));
        offset = offset + current_size;
    end
end

end
