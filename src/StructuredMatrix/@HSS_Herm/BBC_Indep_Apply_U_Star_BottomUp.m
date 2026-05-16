function Y = BBC_Indep_Apply_U_Star_BottomUp(A)

arguments (Input)
    A HSS_Herm;
end

arguments(Output)
    Y (:, :) double;
end

if A.leaf_ == 0
    Y = 0;
    offset = 0;
    for i = 1 : A.num_children_
        child_Y = A.children_{i}.BBC_Indep_Apply_U_Star_BottomUp();
        current_size = size(child_Y, 1);
        Y = Y + A.Umat_((offset + 1) : (offset + current_size), :)' * child_Y;
        offset = offset + current_size;
    end
else
    Y = A.Umat_' * A.BBC_Y_;
    A.BBC_Y_ = [];
end

end