function Apply_Downward(A, level)
% Apply_Downward

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        A.fvec_ = A.Amat_ * A.uvec_ +  A.Umat_ * A.fhvec_;
    else
        for i = 1 : A.num_children_
            A.children_{i}.fhvec_ = A.Rmat_{i} * A.fhvec_;
        end

        for i = 1 : A.num_children_
            for j = 1 : (i - 1)
                A.children_{i}.fhvec_ = A.children_{i}.fhvec_ ...
                    + A.Bmat_{i, j} * A.children_{j}.uhvec_;
            end

            for j = (i + 1) : A.num_children_
                A.children_{i}.fhvec_ = A.children_{i}.fhvec_ ...
                    + A.Bmat_{j, i}' * A.children_{j}.uhvec_;
            end
        end
    end

    % Clear.
    A.uvec_ = [];
    A.uhvec_ = [];
    A.fhvec_ = [];
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.Apply_Downward(level);
    end
end

end