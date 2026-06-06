function OutCompression(BF, tol)

arguments (Input)
    BF Butterfly;
    tol (1, 1) double;
end

L = BF.L_;
h_x = BF.h_x_;
h_xi = BF.h_xi_;
L_x = BF.L_x_;
L_xi = BF.L_xi_;

num_children = BF.num_children_;
ch_list = BF.ch_list_;

% Compress M.
level_x = h_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = h_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});

L_cell = cell(size(BF.M_));
R_cell = cell(size(BF.M_));
for ind_tau = 1 : m_x
    for ind_sigma = 1 : m_xi
        [U, S, V] = MySVDSketch(BF.M_{ind_tau, ind_sigma}, tol);
        L_cell{ind_tau, ind_sigma} = U*sqrt(S);
        BF.M_{ind_tau, ind_sigma} = eye(size(S));
        R_cell{ind_tau, ind_sigma} = sqrt(S)*V';
    end
end

% Compress G.
cnt_G = L_x - h_x + 1;
for level = h_x : (L_x - 1)
    level_x = level + 1;
    ind_level_x = level_x + 1;
    level_x_par = level_x - 1;
    ind_level_x_par = level_x_par + 1;
    m_x = length(BF.tree_{ind_level_x});
    m_x_par = length(BF.tree_{ind_level_x_par});

    level_xi = L - level_x;
    ind_level_xi = level_xi + 1;
    level_xi_ch = level_xi + 1;
    ind_level_xi_ch = level_xi_ch + 1;
    m_xi = length(BF.tree_{ind_level_xi});
    m_xi_ch = length(BF.tree_{ind_level_xi_ch});

    cnt_G = cnt_G - 1;

    % Update G.
    for ind_alpha = 1 : m_x_par
        for ind_sigma = 1 : m_xi
            alpha = BF.tree_{ind_level_x_par}{ind_alpha};

            sigma = BF.tree_{ind_level_xi}{ind_sigma};

            for ch_alpha = ch_list
                ind_tau = num_children * alpha.order_ + ch_alpha + 1;
                tau = BF.tree_{ind_level_x}{ind_tau};

                G_tau_sigma = BF.G_{cnt_G}{ind_tau, ind_sigma};
                G_tau_beta_cell = cell(1, num_children);
                beta_offset = 0;
                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    beta = BF.tree_{ind_level_xi_ch}{ind_beta};

                    L_alpha_beta = L_cell{ind_alpha, ind_beta};

                    beta_size = size(L_alpha_beta, 1);
                    G_tau_beta = G_tau_sigma(:, (beta_offset + 1) : (beta_offset + beta_size));
                    G_tau_beta_cell{ch_sigma + 1} = G_tau_beta * L_alpha_beta;

                    beta_offset = beta_offset + beta_size;
                end

                BF.G_{cnt_G}{ind_tau, ind_sigma} = cell2mat(G_tau_beta_cell);
            end
        end
    end

    % Compress G.
    L_cell = cell(m_x, m_xi);
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            [U, S, V] = MySVDSketch(BF.G_{cnt_G}{ind_tau, ind_sigma}, tol);
            L_cell{ind_tau, ind_sigma} = U * S;
            BF.G_{cnt_G}{ind_tau, ind_sigma} = V';
        end
    end
end

% Update U.
level_x = L_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = L - level_x;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});
for ind_tau = 1 : m_x
    for ind_sigma = 1 : m_xi
        BF.U_{ind_tau, ind_sigma} = BF.U_{ind_tau, ind_sigma} * L_cell{ind_tau, ind_sigma};
    end
end

% Compress H.
cnt_H = 0;
for level = h_xi : (L_xi - 1)
    level_xi = level + 1;
    ind_level_xi = level_xi + 1;
    level_xi_par = level_xi - 1;
    ind_level_xi_par = level_xi_par + 1;
    m_xi = length(BF.tree_{ind_level_xi});
    m_xi_par = length(BF.tree_{ind_level_xi_par});

    level_x = L - level_xi;
    ind_level_x = level_x + 1;
    level_x_ch = level_x + 1;
    ind_level_x_ch = level_x_ch + 1;
    m_x = length(BF.tree_{ind_level_x});
    m_x_ch = length(BF.tree_{ind_level_x_ch});

    cnt_H = cnt_H + 1;

    % Update H.
    for ind_beta = 1 : m_xi_par
        for ind_tau = 1 : m_x
            beta = BF.tree_{ind_level_xi_par}{ind_beta};

            tau = BF.tree_{ind_level_x}{ind_tau};

            for ch_beta = ch_list
                ind_sigma = num_children * beta.order_ + ch_beta + 1;
                sigma = BF.tree_{ind_level_xi}{ind_sigma};

                H_tau_sigma = BF.H_{cnt_H}{ind_tau, ind_sigma};
                H_alpha_sigma_cell = cell(num_children, 1);
                alpha_offset = 0;
                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    alpha = BF.tree_{ind_level_x_ch}{ind_alpha};

                    R_alpha_beta = R_cell{ind_alpha, ind_beta};

                    alpha_size = size(R_alpha_beta, 2);
                    H_alpha_sigma = H_tau_sigma((alpha_offset + 1) : (alpha_offset + alpha_size), :);
                    H_alpha_sigma_cell{ch_tau + 1} = R_alpha_beta * H_alpha_sigma;

                    alpha_offset = alpha_offset + alpha_size;
                end

                BF.H_{cnt_H}{ind_tau, ind_sigma} = cell2mat(H_alpha_sigma_cell);
            end
        end
    end
    
    % Compress H.
    R_cell = cell(m_x, m_xi);
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            [U, S, V] = MySVDSketch(BF.H_{cnt_H}{ind_tau, ind_sigma}, tol);
            R_cell{ind_tau, ind_sigma} = S * V';
            BF.H_{cnt_H}{ind_tau, ind_sigma} = U;
        end
    end
end

% Update V.
level_xi = L_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});

level_x = L - level_xi;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});
for ind_tau = 1 : m_x
    for ind_sigma = 1 : m_xi
        BF.V_{ind_tau, ind_sigma} = R_cell{ind_tau, ind_sigma} * BF.V_{ind_tau, ind_sigma};
    end
end

end