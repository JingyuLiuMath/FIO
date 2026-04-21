function Apply_Root(A)
% Apply_Root

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 1
    A.fvec_ = A.Amat_ * A.uvec_;
else
    for i = 1 : A.num_children_
        A.children_{i}.fhvec_ ...
            = zeros(A.children_{i}.rank_, A.vec_col_size_);
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


end