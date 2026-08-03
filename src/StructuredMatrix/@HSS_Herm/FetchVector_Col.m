function u = FetchVector_Col(A)
% FetchVector_Col

arguments (Input)
    A HSS_Herm;
end

arguments (Output)
    u (:, :) double;
end

u = zeros(A.size_, A.vec_col_size_);
A.vec_col_size_ = [];

if A.leaf_ == 1
    u = A.uvec_;

    % Clear.
    A.uvec_ = [];
else
    % Recursion.
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.size_;
        u((offset + 1) : (offset + current_size), :) ...
            = A.children_{i}.FetchVector_Col();
        offset = offset + current_size;
    end
end

end
