function out = BlackBoxConstruct_Indep(A, op_A, rank_func, tol, verbose)
% BlackBoxConstruct

arguments (Input)
    A HSS_Herm;
    op_A function_handle;
    rank_func function_handle;
    tol (1, 1) double;
    verbose (1, 1) double = 1;
end

arguments (Output)
    out struct;
end

out.s_list = zeros(A.max_level_ + 1, 1);
out.r_list = zeros(A.max_level_ + 1, 1);

out.t_sampling_list = zeros(A.max_level_ + 1, 1);
out.t_construct_generators_list = zeros(A.max_level_ + 1, 1);
out.t_level_list = zeros(A.max_level_ + 1, 1);

p = 5;
s_total = 0;
% Recursive construction.
for level = A.max_level_ : -1 : 0
    if verbose == 1
        fprintf("    \n");
        fprintf("    level: %d\n", level);
    end
    % Settings.
    max_level_size = A.BBC_Indep_LevelSize(level);
    target_rank = ceil(rank_func(level));
    if level ~= 0
        s = target_rank + max_level_size + p;
    else
        s = max_level_size + p;
    end
    s_total = s_total + s;

    t_level_start = tic;
    % Sampling.
    t_sampling_start = tic;
    if level == A.max_level_
        Y = randn(A.global_size_, s);
        A.BBC_FillAuxiliaryMatrix(Y, op_A(Y));
    else
        A.BBC_Indep_TopDown(level + 1, s);
        % A.BBC_Indep_FillY_Leaf(op_A(A.BBC_Indep_FetchY_Leaf(s)));
        Y = A.BBC_Indep_FetchY_Leaf(s);
        Y = op_A(Y);
        A.BBC_Indep_FillY_Leaf(Y);
        A.BBC_Indep_BottomUp(level + 1);
    end
    t_sampling = toc(t_sampling_start);

    t_construct_generators_start = tic;
    if level ~= 0
        level_rank = A.BBC_Indep_ConstructGenerators(...
            level, s, target_rank, tol);
    else
        A.BBC_Indep_ConstructRootGenerators(s);
    end
    t_construct_generators = toc(t_construct_generators_start);

    t_level = toc(t_level_start);

    out.s_list(level + 1) = s;
    if level ~= 0 
        out.r_list(level + 1) = level_rank;
    end
    out.t_sampling_list(level + 1) = t_sampling;
    out.t_construct_generators_list(level + 1) = t_construct_generators;
    out.t_level_list(level + 1) = t_level;

    if verbose == 1
        fprintf("    t_level: %.1e\n", t_level);
        fprintf("    t_sampling: %.1e\n", t_sampling);
        fprintf("    t_construct_generators: %.1e\n", t_construct_generators);
        fprintf("    max_level_size: %d\n", max_level_size);
        if level ~= 0
            fprintf("    target_rank: %d\n", target_rank);
            fprintf("    level_rank: %d\n", level_rank);
        end
        fprintf("    num of samples: %d\n", s);
    end
end

out.s_total = s_total;
fprintf("  total num of samples: %d\n", s_total);


out.t_sampling = sum(out.t_sampling_list);
out.t_construct_generators = sum(out.t_construct_generators_list);

A.BBC_Indep_AssignR();

% Eliminate B{i, i} matrices.
t_postprocessing_start = tic;
A.BBC_EliminateRootBMat();

for level = 1 : 1 : A.max_level_
    A.BBC_EliminateBMat(level);
end
t_postprocessing = toc(t_postprocessing_start);

out.t_postprocessing = t_postprocessing;

if verbose == 1
    fprintf("   t_postprocessing: %.1e\n", t_postprocessing);
end

end