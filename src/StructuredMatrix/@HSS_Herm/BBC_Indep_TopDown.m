function BBC_Indep_TopDown(A, level, s)

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
    s (1, 1) double;
end

if A.level_ == level
    Omega = randn(A.level_size_, s);
    A.BBC_Omega_ = Omega;
    A.BBC_Indep_Apply_U_TopDown(Omega);
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.BBC_Indep_TopDown(level, s);
    end
end

end