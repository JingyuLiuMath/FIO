function BlackBoxConstruct_Indep_New(A, op_A, rank_func, tol, verbose)
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
    if level == A.max_level_
        Y = op_A(Omega);
    else
        A.BBC_Indep_FillY(level + 1, Omega);
        A.BBC_Indep_Apply_U(level + 1);
        Y = A.BBC_Indep_FetchY_Leaf(s);
        Y = op_A(Y);
        A.BBC_Indep_FillY_Leaf(Y);
        Y = [];
        A.BBC_Indep_Apply_U_Star(level + 1);
        Y = A.BBC_Indep_FetchY(level + 1, s);
    end

    if level ~= 0
        A.BBC_Indep_ConstructGenerators_New(...
            level, ...
            Omega, Y, ...
            target_rank, ...
            tol);
    else
        A.BBC_Indep_ConstructRootGenerators_New(Omega, Y);
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
        mem = byte_to_gb(A.Storage());
        mem = mem + byte_to_gb(byte_size(Omega));
        mem = mem + byte_to_gb(byte_size(Y));
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