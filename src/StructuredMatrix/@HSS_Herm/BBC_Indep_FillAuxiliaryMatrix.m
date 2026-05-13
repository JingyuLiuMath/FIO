function BBC_Indep_FillAuxiliaryMatrix(A, level, Omega, Y)
% BBC_FillAuxiliaryMatrix

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
    Omega (:, :) double;
    Y (:, :) double;
end

if A.level_ == level
    A.BBC_Omega_ = Omega;
    A.BBC_Y_ = Y;
else
    % Recursion.
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.level_size_;
        A.children_{i}.BBC_Indep_FillAuxiliaryMatrix(...
            level, ...
            Omega((offset + 1) : (offset + current_size), :), ...
            Y((offset + 1) : (offset + current_size), :));
        offset = offset + current_size;
    end
end

end