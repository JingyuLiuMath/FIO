function BlackBoxConstruct_Indep(A, op_A, rank_func, tol, verbose)
% BlackBoxConstruct

arguments (Input)
    A HSS_Herm;
    op_A function_handle;
    rank_func function_handle;
    tol (1, 1) double;
    verbose (1, 1) double = 1;
end

matlab_mem0 = getMemoryInfo();
fprintf("    initial used memory: %.1e GB\n", matlab_mem0);

p = 5;
s_total = 0;
% Recursive construction.
for level = A.max_level_ : -1 : 0
    % Settings.
    max_level_size = A.BBC_Indep_LevelSize(level);
    target_rank = ceil(rank_func(level));
    if level ~= 0
        s = target_rank + max_level_size + p;
    else
        s = max_level_size + p;
    end
    s_total = s_total + s;

    % Sampling.
    if level == A.max_level_
        Omega = randn(A.global_size_, s);
        A.BBC_FillAuxiliaryMatrix(Omega, op_A(Omega));
        clear Omega;
    else
        A.BBC_Indep_TopDown(level + 1, s);
        A.BBC_Indep_FillY_Leaf(op_A(A.BBC_Indep_FetchY_Leaf(s)));
        A.BBC_Indep_BottomUp(level + 1);
    end

    if level ~= 0
        A.BBC_Indep_ConstructGenerators(...
            level, s, target_rank, tol);
    else
        A.BBC_Indep_ConstructRootGenerators(s);
    end

    if verbose == 1
        fprintf("    \n");
        fprintf("    level: %d\n", level);
        fprintf("    max_level_size: %d\n", max_level_size);
        fprintf("    target_rank: %d\n", target_rank);
        fprintf("    num of samples: %d\n", s);
        mem = byte_to_gb(A.Storage());
        fprintf("    used memory: %.1e GB\n", mem);
        matlab_mem = getMemoryInfo();
        fprintf("    incr used memory: %.1e GB\n", matlab_mem - matlab_mem0);
        fprintf("    total used memory: %.1e GB\n", matlab_mem);
    end
end

fprintf("  total num of samples: %d\n", s_total);

A.BBC_Indep_AssignR();

% Eliminate B{i, i} matrices.
A.BBC_EliminateRootBMat();

for level = 1 : 1 : A.max_level_
    A.BBC_EliminateBMat(level);
end

end