function f = Apply(A, u)
% Apply f = A * u + f.

arguments (Input)
    A HSS_Herm;
    u (:, :) double;
end

arguments (Output)
    f (:, :) double;
end

A.FillVector_Col(u);

% Upward.
for level = A.max_level_ : -1 : 1
    A.Apply_Upward(level);
end

% Root.
A.Apply_Root();

% Downward
for level = 1 : 1 : A.max_level_
    A.Apply_Downward(level);
end

f = A.FetchVector_Row();

end