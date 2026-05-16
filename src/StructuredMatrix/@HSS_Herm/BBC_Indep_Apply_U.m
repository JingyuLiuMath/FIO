function BBC_Indep_Apply_U(A, level)

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
end

if A.level_ == level
    A.BBC_Indep_Apply_U_Downward();
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.BBC_Indep_Apply_U(level);
    end
end

end