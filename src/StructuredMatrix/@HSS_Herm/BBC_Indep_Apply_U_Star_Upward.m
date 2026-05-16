function BBC_Indep_Apply_U_Star_Upward(A)

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    A.BBC_Y_ = [];
    for i = 1 : A.num_children_
        A.children_{i}.BBC_Indep_Apply_U_Star_Upward();
        A.BBC_Y_ = [A.BBC_Y_; A.children_{i}.BBC_Y_];
        A.children_{i}.BBC_Y_ = [];
    end
end

A.BBC_Y_ = A.Umat_' * A.BBC_Y_;

end