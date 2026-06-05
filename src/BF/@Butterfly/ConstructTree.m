function ConstructTree(BF, N)

arguments (Input)
    BF Butterfly;
    N (1, 1) double;
end

root_node = BFNode(N, 0, N - 1, 0, 0);

BF.L_ = ceil(log2(N));
BF.h_x_ = ceil(BF.L_ / 2);
BF.h_xi_ = BF.L_ - BF.h_x_;
fprintf("  L: %d, h_x: %d, h_xi: %d\n", BF.L_, BF.h_x_, BF.h_xi_);

BF.tree_ = cell(1, BF.L_ + 1);
for level = 0 : BF.L_
    BF.tree_{level + 1} = cell(1, 2^level);
end

level = 0;
ind_level = level + 1;
BF.tree_{ind_level}{1} = root_node;
for level = 0 : (BF.L_ - 1)
    ind_level = level + 1;
    for k = 1 : 2^level
        curr_node = BF.tree_{ind_level}{k};
        c1_size = floor(curr_node.ind_size_ / 2);
        c2_size = curr_node.ind_size_ - c1_size;
        ch_size_list = [c1_size, c2_size];
        ch_index_list = [2 * k - 1, 2 * k];
        offset = 0;
        for ck = [0, 1]
            ch_index = ch_index_list(ck + 1);
            ch_size = ch_size_list(ck + 1);
            child_node = BFNode(N, ...
                curr_node.ind_start_ + offset, ...
                curr_node.ind_start_ + offset + ch_size - 1, ...
                level + 1, 2 * curr_node.order_ + ck);
            BF.tree_{ind_level + 1}{ch_index} = child_node;

            offset = offset + ch_size;
        end
    end
end

end