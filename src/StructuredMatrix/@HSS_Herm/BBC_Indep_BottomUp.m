function BBC_Indep_BottomUp(A, level)

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
end

if A.level_ == level
    A.BBC_Y_ = A.BBC_Indep_Apply_U_Star_BottomUp();
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.BBC_Indep_BottomUp(level);
    end
end

end