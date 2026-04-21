function ULV_Solve_Root(A)
% ULV_Solve_Root

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    A.ULV_Solve_MergeVector();
end

A.uvec_ = A.ULV_A_re_re_' \ (A.ULV_A_re_re_ \ A.fvec_);

% Clear.
A.fvec_ = [];

if A.leaf_ == 0
    A.ULV_Solve_SplitVector();
end

end