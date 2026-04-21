function ULV_Factor(A)
% ULV

arguments (Input)
    A HSS_Herm;
end

% Elimination and merge.
for level = A.max_level_ : -1 : 1
    A.ULV_Eliminate(level);
end

A.ULV_RootFactor();

end