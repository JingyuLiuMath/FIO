function u = ULV_Solve(A, f)
% ULV_Solve

arguments (Input)
    A HSS_Herm;
    f (:, :) double;
end

arguments (Output)
    u (:, :) double;
end

A.FillVector_Row(f);

% Upward.
for level = A.max_level_ : -1 : 1
    A.ULV_Solve_Upward(level);
end

% Root.
A.ULV_Solve_Root();

% Downward
for level = 1 : 1 : A.max_level_
    A.ULV_Solve_Downward(level);
end

u = A.FetchVector_Col();

end