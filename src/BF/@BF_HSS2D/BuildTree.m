function BuildTree(G, min_points)

arguments (Input)
    G BF_HSS2D;
    min_points (1, 1) double;
end

if G.size_ <= min_points
    G.leaf_ = 1;
    G.max_level_ = G.level_;

    x_ind = (G.x_freq_start_ : G.x_freq_end_)' + G.nx_ / 2 + 1;
    y_ind = (G.y_freq_start_ : G.y_freq_end_)' + G.ny_ / 2 + 1;
    xy_ind = TensorProduct2D(x_ind, y_ind);
    G.perm_ = sub2ind([G.nx_, G.ny_], xy_ind(:, 1), xy_ind(:, 2));
else
    % Partition.
    G.num_children_ = 4;
    G.children_ = cell(1, G.num_children_);
    x_c1_size = floor(G.x_size_ / 2);
    x_child_size = [x_c1_size, G.x_size_ - x_c1_size];
    y_c1_size = floor(G.y_size_ / 2);
    y_child_size = [y_c1_size, G.y_size_ - y_c1_size];

    it = 0;
    offset = 0;
    y_offset = 0;
    for ity = 1 : 2
        curr_y_size = y_child_size(ity);
        x_offset = 0;
        for itx = 1 : 2
            it = it + 1;
            curr_x_size = x_child_size(itx);
            curr_size = curr_x_size * curr_y_size;
            G.children_{it} = BF_HSS2D(G.nx_, G.ny_, ...
                G.x_freq_start_ + x_offset, ...
                G.x_freq_start_ + x_offset + curr_x_size - 1, ...
                G.y_freq_start_ + y_offset, ...
                G.y_freq_start_ + y_offset + curr_y_size - 1, ...
                G.level_ + 1, ...
                G.offset_ + offset);
            offset = offset + curr_size;
            x_offset = x_offset + curr_x_size;
        end
        y_offset = y_offset + curr_y_size;
    end

    G.perm_ = zeros(G.size_, 1);
    % Recursion.
    offset = 0;
    for i = 1 : G.num_children_
        curr_size = G.children_{i}.size_;
        G.children_{i}.BuildTree(min_points);
        G.perm_((offset + 1) : (offset + curr_size)) = G.children_{i}.perm_;
        G.max_level_ = max(G.max_level_, G.children_{i}.max_level_);
        offset = offset + curr_size;
    end
end

if G.level_ == 0
    [~, G.perm_inv_] = sort(G.perm_, "ascend");
end

end