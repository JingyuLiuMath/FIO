function f = FetchVector_Row(A)
% FetchVector_Row

% Jingyu Liu, November 18, 2024.

arguments (Input)
    A HSS_Herm;
end

arguments (Output)
    f (:, :) double;
end

f = zeros(A.size_, A.vec_col_size_);
A.vec_col_size_ = [];

if A.leaf_ == 1
    f = A.fvec_;

    % Clear.
    A.fvec_ = [];
else
    % Recursion.
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.size_;
        f((offset + 1) : (offset + current_size), :) ...
            = A.children_{i}.FetchVector_Row();
        offset = offset + current_size;
    end
end

end