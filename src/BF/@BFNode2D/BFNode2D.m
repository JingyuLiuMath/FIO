classdef BFNode2D

    properties
        n_ (1, 1) double;
        half_n_ (1, 1) double;
        N_ (1, 1) double;

        node1_ BFNode;
        node2_ BFNode;

        % ind1_start_ (1, 1) double;
        % ind1_end_ (1, 1) double;
        % ind1_size_ (1, 1) double;

        % ind2_start_ (1, 1) double;
        % ind2_end_ (1, 1) double;
        % ind2_size_ (1, 1) double;

        ind_size_ (1, 1) double;

        level_ (1, 1) double;
        order_ (1, 1) double;
    end

    methods
        function B = BFNode2D(n, ...
                ind1_start, ind1_end, ...
                ind2_start, ind2_end, ...
                level, order)
            B.n_ = n;
            B.half_n_ = n / 2;
            B.N_ = n^2;

            % B.ind1_start_ = ind1_start;
            % B.ind1_end_ = ind1_end;
            % B.ind1_size_ = ind1_end - ind1_start + 1;
            B.node1_ = BFNode(n, ind1_start, ind1_end, level, -1);
            B.node2_ = BFNode(n, ind2_start, ind2_end, level, -1);

            % B.ind2_start_ = ind2_start;
            % B.ind2_end_ = ind2_end;
            % B.ind2_size_ = ind2_end - ind2_start + 1;

            B.ind_size_ = B.node1_.ind_size_ * B.node2_.ind_size_;
            
            B.level_ = level;
            B.order_ = order;
        end

    end
end