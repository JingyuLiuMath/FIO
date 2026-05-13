function BlackBoxConstruct_Indep(A, op_A, rank_func, tol, verbose)
% BlackBoxConstruct

arguments (Input)
    A HSS_Herm;
    op_A function_handle;
    rank_func function_handle;
    tol (1, 1) double;
    verbose (1, 1) double = 0;
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
    s = max(r + max_level_size, 2 * r);
    if verbose == 1
        fprintf("  total_level_size: %d\n", total_level_size);
        fprintf("  num of samples: %d\n", s);
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
    
    A.BBC_Indep_FillAuxiliaryMatrix(level, Omega, Y);
    if level ~= 0
        U_level_list = {};
        U_level_list = A.BBC_Indep_ConstructGenerators(...
            level, ...
            r, ...
            U_level_list, ...
            tol);
        U_list{level} = sparse_blkdiag(U_level_list{:});
    else
        A.BBC_Indep_ConstructRootGenerators();
    end
end

fprintf("  total num of samples: %d\n", s_total);

% Eliminate B{i, i} matrices.
A.BBC_EliminateRootBMat();

for level = 1 : 1 : A.max_level_
    A.BBC_EliminateBMat(level);
end

end