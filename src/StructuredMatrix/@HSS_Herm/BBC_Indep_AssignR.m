function BBC_Indep_AssignR(A)

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    % Assign R and W.
    if A.level_ ~= 0
        offset = 0;
        for i = 1 : A.num_children_
            current_size = A.children_{i}.rank_;
            A.Rmat_{i} = A.Umat_(...
                (offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
        end
        A.Umat_ = [];
    end

    for i = 1 : A.num_children_
        A.children_{i}.BBC_Indep_AssignR();
    end
end

end