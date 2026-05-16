function BlackBoxConstruct_Indep(A, op_A, rank_func, tol, verbose)
% BlackBoxConstruct

arguments (Input)
    A HSS_Herm;
    op_A function_handle;
    rank_func function_handle;
    tol (1, 1) double;
    verbose (1, 1) double = 1;
end

p = 5;
s_total = 0;
U_list = cell(1, A.max_level_);
% Recursive construction.
for level = A.max_level_ : -1 : 0
    % Settings.
    max_level_size = A.BBC_Indep_LevelSize(level);
    total_level_size = A.level_size_;
    target_rank = ceil(rank_func(level));
    r = target_rank + p;
    if level ~= 0
        s = r + max_level_size;
    else
        s = p + max_level_size;
    end
    s_total = s_total + s;

    % Sampling.
    Omega = randn(total_level_size, s);
    Y = Omega;
    for it_level = (level + 1) : 1 : A.max_level_
        Y = U_list{it_level} * Y;
    end
    Y = op_A(Y);
    for it_level = A.max_level_ : -1 : (level + 1)
        Y = U_list{it_level}' * Y;
    end

    if level ~= 0
        U_level_list = {};
        U_level_list = A.BBC_Indep_ConstructGenerators(...
            level, ...
            Omega, Y, ...
            target_rank, ...
            U_level_list, ...
            tol);
        U_list{level} = sparse_blkdiag(U_level_list{:});
        U_level_list = {};
    else
        A.BBC_Indep_ConstructRootGenerators(Omega, Y);
    end

    Omega = [];
    Y = [];

    if verbose == 1
        fprintf("    \n");
        fprintf("    level: %d\n", level);
        fprintf("    max_level_size: %d\n", max_level_size);
        fprintf("    total_level_size: %d\n", total_level_size);
        fprintf("    target_rank: %d\n", target_rank);
        fprintf("    num of samples: %d\n", s);
        mem = double_to_gb(A.Storage());
        for it = (level + 1) : A.max_level_
            mem = mem + double_to_gb(nnz(U_list{it}));
        end
        fprintf("    used memory: %.1e GB\n", mem);
        matlab_mem = memory;
        fprintf("    total used memory: %.1e GB\n", matlab_mem.MemUsedMATLAB / 10^9);
    end
end

fprintf("  total num of samples: %d\n", s_total);

% Eliminate B{i, i} matrices.
A.BBC_EliminateRootBMat();

for level = 1 : 1 : A.max_level_
    A.BBC_EliminateBMat(level);
end

end