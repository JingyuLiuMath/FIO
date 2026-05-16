function Y = BBC_Indep_FetchY(A, level, s)

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
    s (1, 1) double;
end

arguments (Output)
    Y (:, :) double;
end

if A.level_ == level
    Y = A.BBC_Y_;
    A.BBC_Y_ = [];
    return;
end

Y = zeros(A.level_size_, s);

% Recursion.
offset = 0;
for i = 1 : A.num_children_
    current_size = A.children_{i}.level_size_;
    Y((offset + 1) : (offset + current_size), :) ...
        = A.children_{i}.BBC_Indep_FetchY(level, s);
    offset = offset + current_size;
end

end