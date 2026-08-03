function ULV_Solve_Downward(A, level)
% ULV_Solve_Downward

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
end

if A.level_ == level
    % Update f_re.
    A.ULV_u_re_ = A.ULV_u_re_ - A.ULV_A_sk_re_' * A.ULV_u_sk_;
    A.ULV_u_re_ = A.ULV_A_re_re_' \ A.ULV_u_re_;
    A.uvec_ = A.ULV_Q_ * [A.ULV_u_re_; A.ULV_u_sk_];

    % Clear.
    A.ULV_u_re_ = [];
    A.ULV_u_sk_ = [];

    if A.leaf_ == 0
        A.ULV_Solve_SplitVector();
    end
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.ULV_Solve_Downward(level);
    end
end

end
