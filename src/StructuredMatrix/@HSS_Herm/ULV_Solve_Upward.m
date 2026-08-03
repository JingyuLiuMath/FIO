function ULV_Solve_Upward(A, level)
% ULV_Solve_Upward

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 0
        A.ULV_Solve_MergeVector();
    end

    A.fvec_ = A.ULV_Q_' * A.fvec_;
    A.ULV_u_re_ = A.ULV_A_re_re_ \  A.fvec_(1 : A.ULV_re_size_, :);
    A.ULV_f_sk_ = A.fvec_((A.ULV_re_size_ + 1) : end, :) ...
        - A.ULV_A_sk_re_ * A.ULV_u_re_;

    % Clear.
    A.fvec_ = [];
elseif A.leaf_ == 0
    % Recursion.
    for i = 1 : A.num_children_
        A.children_{i}.ULV_Solve_Upward(level);
    end
end

end
