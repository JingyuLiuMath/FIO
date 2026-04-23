function ULV_RootFactor(A)
% ULV_RootFactor

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    A.ULV_Merge();
end

A.Amat_ = (A.Amat_ + A.Amat_') / 2;
[A.ULV_A_re_re_, flag] = chol(A.Amat_, 'lower');
A.ULV_lu_root_ = 0;
if flag ~= 0
    [A.ULV_L_root_, A.ULV_U_root_] = lu(A.Amat_);
    A.ULV_lu_root_ = 1;
end

% Clear data.
A.Amat_ = [];

end