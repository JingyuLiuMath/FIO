function ULV_Eliminate(A, level)
% ULV_Eliminate

arguments (Input)
    A HSS_Herm;
    level (1, 1) double;
end

if A.level_ == level
    if A.leaf_ == 0
        A.ULV_Merge();
    end
    
    m = size(A.Umat_, 1);
    k = A.rank_;
    
    % Step 1. Zeroing out U.
    [A.ULV_Q_, R] = qr(A.Umat_);
    reverse_order = m : -1 : 1;
    A.ULV_Q_ = A.ULV_Q_(:, reverse_order);
    R = R(reverse_order, :);
    sk_size = k;
    re_size = m - k;
    A.ULV_re_size_ = re_size;
    A.ULV_sk_size_ = sk_size;
    A.ULV_U_sk_ = R((re_size + 1) : m, :);
    A.Amat_ = A.ULV_Q_' * A.Amat_ * A.ULV_Q_;
    A.Amat_ = (A.Amat_ + A.Amat_') / 2;

    % Step 2. Cholesky on A.
    % chol only uses the triangular part of the matrix.
    A.ULV_A_re_re_ = chol(A.Amat_(1 : re_size, 1 : re_size), 'lower');
    A.ULV_A_sk_re_ = A.Amat_((1 + re_size) : end, 1 : re_size) / A.ULV_A_re_re_';
    A.ULV_A_sk_sk_ = A.Amat_((1 + re_size) : end, (1 + re_size) : end) ...
        - A.ULV_A_sk_re_ * A.ULV_A_sk_re_';

    % Clear data.
    A.Amat_ = [];
    A.Umat_ = [];
elseif A.leaf_ == 0
    for i = 1 : A.num_children_
        A.children_{i}.ULV_Eliminate(level);
    end
end

end