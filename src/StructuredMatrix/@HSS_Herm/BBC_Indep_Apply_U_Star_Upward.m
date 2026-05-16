function BBC_Indep_Apply_U_Star_Upward(A)

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    tmp = zeros(A.rank_, 0);
    offset = 0;
    for i = 1 : A.num_children_
        A.children_{i}.BBC_Indep_Apply_U_Star_Upward();
        child_Y = A.children_{i}.BBC_Y_;
        A.children_{i}.BBC_Y_ = [];
        if isempty(tmp)
            tmp = zeros(A.rank_, size(child_Y, 2));
        end
        current_size = size(child_Y, 1);
        tmp = tmp + A.Umat_((offset + 1) : (offset + current_size), :)' * child_Y;
        offset = offset + current_size;
    end
    A.BBC_Y_ = tmp;
else
    A.BBC_Y_ = A.Umat_' * A.BBC_Y_;
end

end