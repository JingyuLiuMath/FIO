function Apply_Upward(A, level)
% Apply_Upward

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        A.uhvec_ = A.Umat_' * A.uvec_;
    else
        A.uhvec_ = zeros(A.rank_, A.vec_col_size_);
        for i = 1 : A.num_children_
            A.uhvec_ = A.uhvec_ + A.Rmat_{i}' * A.children_{i}.uhvec_;
        end
    end
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.Apply_Upward(level);
    end
end

end