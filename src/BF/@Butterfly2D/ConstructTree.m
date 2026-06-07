function ConstructTree(BF, n_leaf)

arguments (Input)
    BF Butterfly2D;
    n_leaf (1, 1) double;
end


n = BF.n_;
N = BF.N_;
BF.L_ = ceil(log2(n));
BF.h_x_ = ceil(BF.L_ / 2);
BF.h_xi_ = BF.L_ - BF.h_x_;
BF.L_x_ = max(BF.h_x_, ceil(log2(n / n_leaf)));
BF.L_xi_ = max(BF.h_xi_, ceil(log2(n / n_leaf)));
fprintf("  L: %d, h_x: %d, h_xi: %d, L_x: %d, L_xi: %d\n", ...
    BF.L_, BF.h_x_, BF.h_xi_, BF.L_x_, BF.L_xi_);

num_children = 4;
ch_list = [0, 1, 2, 3];
ch_dir_list = [0, 1];

BF.num_children_ = num_children;
BF.ch_list_ = ch_list;

L_max = max(BF.L_x_, BF.L_xi_);

root_node = BFNode2D(n, 0, n - 1, 0, n - 1, 0, 0);

BF.tree_ = cell(1, L_max + 1);
for level = 0 : L_max
    BF.tree_{level + 1} = cell(1, num_children^level);
end

level = 0;
ind_level = level + 1;
BF.tree_{ind_level}{1} = root_node;
for level = 0 : (L_max - 1)
    ind_level = level + 1;
    for k = 1 : num_children^level
        curr_node = BF.tree_{ind_level}{k};

        c11_size = floor(curr_node.node1_.ind_size_ / 2);
        c12_size = curr_node.node1_.ind_size_ - c11_size;
        ch1_size_list = [c11_size, c12_size];

        c21_size = floor(curr_node.node2_.ind_size_ / 2);
        c22_size = curr_node.node2_.ind_size_ - c21_size;
        ch2_size_list = [c21_size, c22_size];

        ch_index_list = num_children * k - flip(ch_list);
        ck = 0;
        c2_offset = 0;
        for ck2 = ch_dir_list
            ch2_size = ch2_size_list(ck2 + 1);
            c1_offset = 0;
            for ck1 = ch_dir_list
                ch_index = ch_index_list(ck + 1);
                ch1_size = ch1_size_list(ck1 + 1);
                child_node = BFNode2D(n, ...
                    curr_node.node1_.ind_start_ + c1_offset, ...
                    curr_node.node1_.ind_start_ + c1_offset + ch1_size - 1, ...
                    curr_node.node2_.ind_start_ + c2_offset, ...
                    curr_node.node2_.ind_start_ + c2_offset + ch2_size - 1, ...
                    level + 1,num_children * curr_node.order_ + ck);
                BF.tree_{ind_level + 1}{ch_index} = child_node;
                ck = ck + 1;

                c1_offset = c1_offset + ch1_size;
            end
            c2_offset = c2_offset + ch2_size;
        end
    end
end

% Assign perm.
BF.x_perm_ = zeros(N, 1);
ind_level_x = BF.L_x_ + 1;
x_offset = 0;
for k = 1 : num_children^BF.L_x_
    ind = BF.tree_{ind_level_x}{k}.Indices();
    x_size = size(ind, 1);
    BF.x_perm_((x_offset + 1) : (x_offset + x_size)) = ind;
    x_offset = x_offset + x_size;
end
[~, BF.x_perm_inv_] = sort(BF.x_perm_, "ascend");


BF.xi_perm_ = zeros(N, 1);
ind_level_xi = BF.L_xi_ + 1;
xi_offset = 0;
for k = 1 : num_children^BF.L_xi_
    ind = BF.tree_{ind_level_xi}{k}.Indices();
    xi_size = size(ind, 1);
    BF.xi_perm_((xi_offset + 1) : (xi_offset + xi_size)) = ind;
    xi_offset = xi_offset + xi_size;
end
[~, BF.xi_perm_inv_] = sort(BF.xi_perm_, "ascend");

end