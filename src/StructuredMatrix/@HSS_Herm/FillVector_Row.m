function FillVector_Row(A, f)
% FillVector_Row

arguments (Input)
    A HSS_Herm;
    f (:, :) double;
end

A.vec_col_size_ = size(f, 2);

if A.leaf_ == 1
    A.fvec_ = f;
else
    % Recursion.
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.size_;
        A.children_{i}.FillVector_Row(...
            f((offset + 1) : (offset + current_size), :));
        offset = offset + current_size;
    end
end

end
