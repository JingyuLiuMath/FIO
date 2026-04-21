function ULV_Merge(A)
% ULV_Merge


arguments (Input)
    A HSS_Herm;
end

sk_size = 0;
for i = 1 : A.num_children_
    sk_size = sk_size + A.children_{i}.ULV_sk_size_;
end

% Merge A.
A.Amat_ = zeros(sk_size, sk_size);
row_offset = 0;
for i = 1 : A.num_children_
    current_row_size = A.children_{i}.ULV_sk_size_;
    col_offset = 0;

    for j = 1 : (i - 1)
        current_col_size = A.children_{j}.ULV_sk_size_;
        A.Amat_((row_offset + 1) : (row_offset + current_row_size), ...
            (col_offset + 1) : (col_offset + current_col_size)) ...
            = A.children_{i}.ULV_U_sk_ ...
            * A.Bmat_{i, j} ...
            * A.children_{j}.ULV_U_sk_';
        col_offset = col_offset + current_col_size;
    end

    j = i;
    current_col_size = A.children_{j}.ULV_sk_size_;
    A.Amat_((row_offset + 1) : (row_offset + current_row_size), ...
        (col_offset + 1) : (col_offset + current_col_size)) ...
        = A.children_{i}.ULV_A_sk_sk_;
    col_offset = col_offset + current_col_size;

    for j = (i + 1) : A.num_children_
        current_col_size = A.children_{j}.ULV_sk_size_;
        A.Amat_((row_offset + 1) : (row_offset + current_row_size), ...
            (col_offset + 1) : (col_offset + current_col_size)) ...
            = A.children_{i}.ULV_U_sk_ ...
            * A.Bmat_{j, i}' ...
            * A.children_{j}.ULV_U_sk_';
        col_offset = col_offset + current_col_size;
    end

    row_offset = row_offset + current_row_size;
end

if A.level_ ~= 0
    % Merge U.
    A.Umat_ = zeros(sk_size, A.rank_);
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.ULV_sk_size_;
        A.Umat_((offset + 1) : (offset + current_size), :) ...
            = A.children_{i}.ULV_U_sk_ * A.Rmat_{i};
        offset = offset + current_size;
    end
end

end