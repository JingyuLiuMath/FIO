function BuildTree(G, min_points)

arguments (Input)
    G BF_HSS;
    min_points (1, 1) double;
end

if G.size_ <= min_points
    G.leaf_ = 1;
    G.max_level_ = G.level_;
else
    % Partition.
    G.num_children_ = 2;
    G.children_ = cell(1, G.num_children_);
    c1_size = floor(G.size_ / 2);
    child_size = [c1_size, G.size_ - c1_size];

    offset = 0;
    for i = 1 : G.num_children_
        current_size = child_size(i);
        G.children_{i} = BF_HSS(...
            G.global_size_, ...
            G.freq_start_ + offset, ...
            G.freq_start_ + offset + current_size - 1, ...
            G.level_ + 1, ...
            G.offset_ + offset);
        offset = offset + current_size;
    end

    % Recursion.
    for i = 1 : G.num_children_
        G.children_{i}.BuildTree(min_points);
        G.max_level_ = max(G.max_level_, G.children_{i}.max_level_);
    end
end

end