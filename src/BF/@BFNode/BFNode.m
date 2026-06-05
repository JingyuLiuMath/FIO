classdef BFNode

    properties
        N_ (1, 1) double;
        half_N_ (1, 1) double;
        ind_start_ (1, 1) double;
        ind_end_ (1, 1) double;
        ind_size_ (1, 1) double;

        level_ (1, 1) double;
        order_ (1, 1) double;
    end

    methods
        function B = BFNode(N, ind_start, ind_end, level, order)
            B.N_ = N;
            B.half_N_ = N / 2;
            B.ind_start_ = ind_start;
            B.ind_end_ = ind_end;
            B.ind_size_ = ind_end - ind_start + 1;
            B.level_ = level;
            B.order_ = order;
        end

    end
end