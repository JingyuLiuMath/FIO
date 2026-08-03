function BlackBoxConstruct(A, op_A, rank_func, tol, verbose)
% BlackBoxConstruct

arguments (Input)
    A HSS_Herm;
    op_A function_handle;
    rank_func function_handle;
    tol (1, 1) double;
    verbose (1, 1) double = 1;
end

% Settings.
leaf_size = A.MaxLeafSize();
mode = 2;
if mode == 1
    total_target_rank = rank_func(1);
else
    total_target_rank = 0;
    for level = 1 : A.max_level_
        total_target_rank = total_target_rank + rank_func(level);
    end
end
p = 5;
s = total_target_rank + p + leaf_size;
fprintf("  sampling mode: %d\n", mode);
fprintf("  total num of samples: %d\n", s);

% Sampling.
Omega = randn(A.size_, s);
Y = op_A(Omega);

A.BBC_FillAuxiliaryMatrix(Omega, Y);

% Recursive construction.
for level = A.max_level_ : -1 : 1
    if verbose == 1
        fprintf("    \n");
        fprintf("    level: %d\n", level);
    end

    target_rank = ceil(rank_func(level));
    t_level_start = tic;
    level_rank = A.BBC_ConstructGenerators(level, target_rank, tol);
    t_level = toc(t_level_start);

    if verbose == 1
        fprintf("    t_level: %.1e\n", t_level);
        fprintf("    target_rank: %d\n", target_rank);
        fprintf("    level_rank: %d\n", level_rank);
    end
end

A.BBC_ConstructRootGenerators();

% Eliminate B{i, i} matrices.
A.BBC_EliminateRootBMat();

for level = 1 : 1 : A.max_level_
    A.BBC_EliminateBMat(level);
end

end
