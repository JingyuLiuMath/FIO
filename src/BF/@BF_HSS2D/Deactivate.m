function inact = Deactivate(A, level, inact)

arguments (Input)
    A BF_HSS2D;
    level (1, 1) double;
    inact (:, 1) double;
end

arguments (Output)
    inact (:, 1) double;
end

if A.level_ == level
    inact = [inact; A.re_];
    A.re_ = [];
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        inact = A.children_{i}.Deactivate(level, inact);
    end
end

end