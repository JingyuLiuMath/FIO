function ULV_Solve_Root(A)
% ULV_Solve_Root

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    A.ULV_Solve_MergeVector();
end

if A.ULV_lu_root_ == 0
    A.uvec_ = A.ULV_A_re_re_' \ (A.ULV_A_re_re_ \ A.fvec_);
else
    A.uvec_ = A.ULV_U_root_ \ (A.ULV_L_root_ \ A.fvec_);
end

% Clear.
A.fvec_ = [];

if A.leaf_ == 0
    A.ULV_Solve_SplitVector();
end

end