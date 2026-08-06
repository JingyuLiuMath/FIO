function BF_adjoint = Adjoint(BF, alpha)
% Adjoint constructs the exact BHP representation of a scaled adjoint.

arguments (Input)
    BF Butterfly2D;
    alpha (1, 1) double = 1;
end

arguments (Output)
    BF_adjoint Butterfly2D;
end

assert(BF.constructed_, "Butterfly2D:Adjoint:NotConstructed", ...
    "Call Construct or BlackBoxConstruct before Adjoint.");

BF_adjoint = Butterfly2D(BF.n_, BF.n_leaf_, 0);

BF_adjoint.h_x_ = BF.h_xi_;
BF_adjoint.h_xi_ = BF.h_x_;
BF_adjoint.L_x_ = BF.L_xi_;
BF_adjoint.L_xi_ = BF.L_x_;

BF_adjoint.x_perm_ = BF.xi_perm_;
BF_adjoint.x_perm_inv_ = BF.xi_perm_inv_;
BF_adjoint.xi_perm_ = BF.x_perm_;
BF_adjoint.xi_perm_inv_ = BF.x_perm_inv_;

BF_adjoint.U_ = AdjointBlockCell(BF.V_);
BF_adjoint.V_ = AdjointBlockCell(BF.U_);
BF_adjoint.M_ = AdjointBlockCell(BF.M_);

num_G = length(BF.G_);
BF_adjoint.H_ = cell(1, num_G);
for k = 1 : num_G
    BF_adjoint.H_{num_G - k + 1} = ...
        AdjointBlockCell(BF.G_{k});
end

num_H = length(BF.H_);
BF_adjoint.G_ = cell(1, num_H);
for k = 1 : num_H
    BF_adjoint.G_{num_H - k + 1} = ...
        AdjointBlockCell(BF.H_{k});
end

BF_adjoint.constructed_ = true;
for i = 1 : numel(BF_adjoint.M_)
    BF_adjoint.M_{i} = alpha * BF_adjoint.M_{i};
end

end

function C_adjoint = AdjointBlockCell(C)

C_adjoint = cell(size(C, 2), size(C, 1));
for i = 1 : size(C, 1)
    for j = 1 : size(C, 2)
        C_adjoint{j, i} = C{i, j}';
    end
end

end
