function ConstructGenerators_ID_Full(G, ...
    exp_phi_func, ...
    x, xi, ...
    inact, ...
    level, rank_or_tol)

arguments (Input)
    G BF_HSS2D;
    exp_phi_func function_handle;
    x (:, 2) double;
    xi (:, 2) double;
    inact (:, 1) double;
    level (1, 1) double;
    rank_or_tol (1, 1) double;
end

N = G.global_size_;
nx = G.nx_;
ny = G.ny_;

if G.level_ == level
    if G.leaf_ == 1
        G.ind_ = G.offset_ + (1 : G.size_)';
    else
        % Merge data from children and construct B.
        G.ind_ = [];
        for i = 1 : G.num_children_
            G.ind_ = [G.ind_; G.children_{i}.ind_];
            c_xi_I = xi(G.children_{i}.ind_, :);
            for j = [1 : (i - 1), (i + 1) : G.num_children_]
                c_xi_J = xi(G.children_{j}.ind_, :);
                G.Bmat_{i, j} = exp_phi_func(x, c_xi_I)' * exp_phi_func(x, c_xi_J);
            end
        end
    end

    % Construct U.
    xi_I = xi(G.ind_, :);
    Ic = (1 : N)';
    self_ind = G.offset_ + (1 : G.size_)';
    Ic([self_ind; inact]) = [];
    xi_Ic = xi(Ic, :);
    A_I_Ic = exp_phi_func(x, xi_I)' * exp_phi_func(x, xi_Ic);
    [sk, U, G.rank_, re] = LowRank_Row_ID(A_I_Ic, rank_or_tol);

    if G.leaf_ == 1
        % Assign U and V.
        G.Umat_ = U;

        % Assign full mat.
        G.Amat_ = exp_phi_func(x, xi_I)' * exp_phi_func(x, xi_I);
    else
        % Assign R and W.
        offset = 0;
        for i = 1 : G.num_children_
            current_size = size(G.children_{i}.ind_, 1);
            G.Rmat_{i} = U((offset + 1) : (offset + current_size), :);
            offset = offset + current_size;
        end

        % Clear data.
        for i = 1 : G.num_children_
            G.children_{i}.ind_ = [];
        end
    end

    % Update row and col.
    G.re_ = G.ind_(re);
    G.ind_ = G.ind_(sk);
else
    for i = 1 : G.num_children_
        G.children_{i}.ConstructGenerators_ID_Full(...
            exp_phi_func, ...
            x, xi, ...
            inact, ...
            level, rank_or_tol);
    end
end

end