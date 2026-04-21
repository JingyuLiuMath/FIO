function BBC_EliminateRootBMat(A)
% BBC_EliminateRootBMat

% Jingyu Liu, December 6, 2024.

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.BBC_A_ = A.Bmat_{i, i};
        A.Bmat_{i, i} = [];
    end
end

end