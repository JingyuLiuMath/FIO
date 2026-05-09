function ConstructRootGenerators2(G, exp_phi_func, x, xi)

arguments (Input)
    G BF_HSS2D;
    exp_phi_func function_handle;
    x (:, 2) double;
    xi (:, 2) double;
end

N = G.global_size_;
nx = G.nx_;
ny = G.ny_;

if G.leaf_ == 1
    G.ind_ = G.offset_ + (1 : G.size_)';

    xi_I = xi(G.ind_, :);

    % Assign full mat.
    G.Amat_ = exp_phi_func(x, xi_I)' * exp_phi_func(x, xi_I);

    % Clear.
    G.ind_ = [];
else
    % Construct B.
    for i = 1 : G.num_children_
        c_xi_I = xi(G.children_{i}.ind_, :);
        for j = [1 : (i - 1), (i + 1) : G.num_children_]
            c_xi_J = xi(G.children_{j}.ind_, :);
            G.Bmat_{i, j} = exp_phi_func(x, c_xi_I)' * exp_phi_func(x, c_xi_J);
        end
    end

    % Clear data.
    for i = 1 : G.num_children_
        G.children_{i}.ind_ = [];
    end
end

end