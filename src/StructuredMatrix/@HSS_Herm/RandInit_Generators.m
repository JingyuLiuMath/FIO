function RandInit_Generators(A, level, hss_rank_rule)

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
    hss_rank_rule string;
end

switch hss_rank_rule
    case "const"
        r = 10;
    case "log"
        r = ceil(log2(A.size_));
    case "sqrt"
        r = ceil(sqrt(A.size_));
end

if A.level_ == level
    A.rank_ =  r;
    if A.leaf_ == 1
        A.Amat_ = randn(A.size_) + randn(A.size_) * 1i;
        A.Amat_ = (A.Amat_ + A.Amat_') / 2;
        A.Umat_ = randn(A.size_, A.rank_);
    else
        for i = 1 : A.num_children_
            A.Rmat_{i} = ...
                randn(A.children_{i}.rank_, A.rank_) ...
                + randn(A.children_{i}.rank_, A.rank_) * 1i;
        end

        for i = 1 : A.num_children_
            for j = 1 : i
                A.Bmat_{i, j} = ...
                    randn(A.children_{i}.rank_, A.children_{j}.rank_) ...
                    + randn(A.children_{i}.rank_, A.children_{j}.rank_) * 1i;
            end
        end
    end
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.RandInit_Generators(level, hss_rank_rule);
    end
end

end