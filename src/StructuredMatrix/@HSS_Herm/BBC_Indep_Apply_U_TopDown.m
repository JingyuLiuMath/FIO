function BBC_Indep_Apply_U_TopDown(A, Y)

arguments (Input)
    A HSS_Herm;
    Y (:, :) double;
end

if A.leaf_ == 0
    Y = A.Umat_ * Y;
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.rank_;
        A.children_{i}.BBC_Indep_Apply_U_TopDown(...
            Y((offset + 1) : (offset + current_size), :));
        offset = offset + current_size;
    end
else
    A.BBC_Y_ = A.Umat_ * Y;
end

end