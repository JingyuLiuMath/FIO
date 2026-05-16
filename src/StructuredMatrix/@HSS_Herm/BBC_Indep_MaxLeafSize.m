function max_leaf_size = BBC_Indep_MaxLeafSize(A)
% MaxLeafSize

arguments (Input)
    A HSS_Herm;
end

arguments (Output)
    max_leaf_size (1, 1) double;
end

if A.leaf_ == 1
    A.level_size_ = A.size_;
    max_leaf_size = A.size_;
else
    A.level_size_ = A.size_;
    max_leaf_size = 0;
    for i = 1 : A.num_children_
        c_max_leaf_size = A.children_{i}.BBC_Indep_MaxLeafSize();
        max_leaf_size = max(max_leaf_size, c_max_leaf_size);
    end
end

end