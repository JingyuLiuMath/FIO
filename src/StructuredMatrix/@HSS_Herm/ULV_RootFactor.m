function ULV_RootFactor(A)
% ULV_RootFactor

arguments (Input)
    A HSS_Herm;
end

if A.leaf_ == 0
    A.ULV_Merge();
end

A.Amat_ = (A.Amat_ + A.Amat_') / 2;
[A.ULV_A_re_re_, ~] = chol(A.Amat_, 'lower');

% Clear data.
A.Amat_ = [];

end