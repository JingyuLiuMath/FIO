function max_level_size = BBC_Indep_LevelSize(A, level)
% MaxLeafSize

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
end

arguments (Output)
    max_level_size (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 1
        A.level_size_ = A.size_;
    else
        A.level_size_ = 0;
        for i = 1 : A.num_children_
            % A.children_{i}.level_size_ = A.children_{i}.rank_;
            A.level_size_ = A.level_size_ + A.children_{i}.rank_;
        end
    end
    max_level_size = A.level_size_;
else
    max_level_size = 0;
    A.level_size_ = 0;
    for i = 1 : A.num_children_
        c_max_level_size = A.children_{i}.BBC_Indep_LevelSize(level);
        A.level_size_ = A.level_size_ + A.children_{i}.level_size_;
        max_level_size = max(max_level_size, c_max_level_size);
    end
end

end