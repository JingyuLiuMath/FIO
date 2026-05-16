function BBC_Indep_Apply_U_Downward(A)

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    tmp = A.Umat_ * A.BBC_Y_;
    A.BBC_Y_ = [];
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.rank_;
        A.children_{i}.BBC_Y_ = tmp((offset + 1) : (offset + current_size), :);
        A.children_{i}.BBC_Indep_Apply_U_Downward();
        offset = offset + current_size;
    end
else
    A.BBC_Y_ = A.Umat_ * A.BBC_Y_;
end

end