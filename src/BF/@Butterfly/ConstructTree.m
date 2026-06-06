function ConstructTree(BF, N, min_points)

arguments (Input)
    BF Butterfly;
    N (1, 1) double;
    min_points (1, 1) double;
end

root_node = BFNode(N, 0, N - 1, 0, 0);

BF.L_ = ceil(log2(N));
BF.h_x_ = ceil(BF.L_ / 2);
BF.h_xi_ = BF.L_ - BF.h_x_;
BF.L_x_ = max(BF.h_x_, ceil(log2(N / min_points)));
BF.L_xi_ = max(BF.h_xi_, ceil(log2(N / min_points)));
fprintf("  L: %d, h_x: %d, h_xi: %d, L_x: %d, L_xi: %d\n", ...
    BF.L_, BF.h_x_, BF.h_xi_, BF.L_x_, BF.L_xi_);

num_children = 2;
ch_list = [0, 1];

BF.num_children_ = num_children;
BF.ch_list_ = ch_list;

BF.tree_ = cell(1, BF.L_ + 1);
for level = 0 : BF.L_
    BF.tree_{level + 1} = cell(1, num_children^level);
end

level = 0;
ind_level = level + 1;
BF.tree_{ind_level}{1} = root_node;
for level = 0 : (BF.L_ - 1)
    ind_level = level + 1;
    for k = 1 : num_children^level
        curr_node = BF.tree_{ind_level}{k};
        c1_size = floor(curr_node.ind_size_ / 2);
        c2_size = curr_node.ind_size_ - c1_size;
        ch_size_list = [c1_size, c2_size];
        ch_index_list = num_children * k - flip(ch_list);
        offset = 0;
        for ck = ch_list
            ch_index = ch_index_list(ck + 1);
            ch_size = ch_size_list(ck + 1);
            child_node = BFNode(N, ...
                curr_node.ind_start_ + offset, ...
                curr_node.ind_start_ + offset + ch_size - 1, ...
                level + 1, num_children * curr_node.order_ + ck);
            BF.tree_{ind_level + 1}{ch_index} = child_node;

            offset = offset + ch_size;
        end
    end
end

end