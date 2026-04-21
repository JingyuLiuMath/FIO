function BlackBoxConstruct(A, op_A, leaf_size, target_rank, tol)
% BlackBoxConstruct

% Jingyu Liu, December 4, 2024.

arguments (Input)
    A HSS_Herm;
    op_A function_handle;
    leaf_size (1, 1) double;
    target_rank (1, 1) double;
    tol (1, 1) double;
end

% Settings.
p = 5;
r = target_rank + p;
s = max(r + leaf_size, 3 * r);

% Sampling.
Omega = randn(A.size_, s);
Y = op_A(Omega);

A.BBC_FillAuxiliaryMatrix(Omega, Y);

% Recursive construction.
for level = A.max_level_ : -1 : 1
    A.BBC_ConstructGenerators(level, r, tol);
end

A.BBC_ConstructRootGenerators();

% Eliminate B{i, i} matrices.
A.BBC_EliminateRootBMat();

for level = 1 : 1 : A.max_level_
    A.BBC_EliminateBMat(level);
end

end