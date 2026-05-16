function BBC_Indep_FillY(A, level, Y)

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
    Y (:, :) double;
end

if A.level_ == level
    A.BBC_Y_ = Y;
elseif A.leaf_ == 0
    % Recursion.
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.level_size_;
        A.children_{i}.BBC_Indep_FillY(...
            level, ...
            Y((offset + 1) : (offset + current_size), :));
        offset = offset + current_size;
    end
end

end