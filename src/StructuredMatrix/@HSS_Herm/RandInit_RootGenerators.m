function RandInit_RootGenerators(A)

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 1
    A.Amat_ = randn(A.size_) + randn(A.size_) * 1i;
    A.Amat_ = (A.Amat_ + A.Amat_') / 2;
else
    for i = 1 : A.num_children_
        for j = 1 : i
            A.Bmat_{i, j} = ...
                randn(A.children_{i}.rank_, A.children_{j}.rank_) ...
                + randn(A.children_{i}.rank_, A.children_{j}.rank_) * 1i;
        end
    end
end

end