function Y = BBC_Indep_FetchY_Leaf(A, s)

arguments (Input)
    A HSS_Herm;
    s (1, 1) double;
end

arguments (Output)
    Y (:, :) double;
end

if A.leaf_ == 1
    Y = A.BBC_Y_;
    A.BBC_Y_ = [];
    return;
end

Y = zeros(A.size_, s);

% Recursion.
offset = 0;
for i = 1 : A.num_children_
    current_size = A.children_{i}.size_;
    Y((offset + 1) : (offset + current_size), :) ...
        = A.children_{i}.BBC_Indep_FetchY_Leaf(s);
    offset = offset + current_size;
end

end