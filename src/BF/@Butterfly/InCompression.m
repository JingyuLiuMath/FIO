function InCompression(BF, tol)

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

% Compress U.
level_x = L_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = L - level_x;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});
L_cell = cell(m_x, m_xi);
for ind_alpha = 1 : m_x
    for ind_beta = 1 : m_xi
        [U, S, V] = MySVDSketch(BF.U_{ind_alpha, ind_beta}, tol);
        BF.U_{ind_alpha, ind_beta} = U;
        L_cell{ind_alpha, ind_beta} = S * V';
    end
end

% Compress G.
cnt_G = 0;
for level = (L_x - 1) : -1 : h_x
    level_x = level;
    ind_level_x = level_x + 1;
    level_x_ch = level_x + 1;
    ind_level_x_ch = level_x_ch + 1;
    m_x = length(BF.tree_{ind_level_x});
    m_x_ch = length(BF.tree_{ind_level_x_ch});

    level_xi = L - level_x;
    ind_level_xi = level_xi + 1;
    level_xi_par = level_xi - 1;
    ind_level_xi_par = level_xi_par + 1;
    m_xi = length(BF.tree_{ind_level_xi});
    m_xi_par = length(BF.tree_{ind_level_xi_par});

    cnt_G = cnt_G + 1;

    % Update G.
    for ind_alpha = 1 : m_x_ch
        for ind_beta = 1 : m_xi_par
            BF.G_{cnt_G}{ind_alpha, ind_beta} = L_cell{ind_alpha, ind_beta} * BF.G_{cnt_G}{ind_alpha, ind_beta};
        end
    end

    % Compress G.
    L_cell = cell(m_x, m_xi);
    for ind_beta = 1 : m_xi_par
        for ind_tau = 1 : m_x
            tau = BF.tree_{ind_level_x}{ind_tau};

            beta = BF.tree_{ind_level_xi_par}{ind_beta};

            U_alpha_sigma_cell = cell(num_children, num_children);
            sigma_offset = 0;
            for ch_beta = ch_list
                ind_sigma = num_children * beta.order_ + ch_beta + 1;
                sigma = BF.tree_{ind_level_xi}{ind_sigma};

                if cnt_G + 1 > length(BF.G_)
                    sigma_size = size(BF.M_{ind_tau, ind_sigma}, 1);
                else
                    sigma_size = size(BF.G_{cnt_G + 1}{ind_tau, ind_sigma}, 1);
                end

                G_alpha_sigma_cell = cell(num_children, 1);
                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    alpha = BF.tree_{ind_level_x_ch}{ind_alpha};

                    G_alpha_beta = BF.G_{cnt_G}{ind_alpha, ind_beta};
                    G_alpha_sigma = G_alpha_beta(:, (sigma_offset + 1) : (sigma_offset + sigma_size));
                    G_alpha_sigma_cell{ch_tau + 1} = G_alpha_sigma;
                end
                G_tau_sigma = cell2mat(G_alpha_sigma_cell);
                [U, S, V] = MySVDSketch(G_tau_sigma, tol);
                L_cell{ind_tau, ind_sigma} = S * V';

                alpha_offset = 0;
                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    alpha = BF.tree_{ind_level_x_ch}{ind_alpha};

                    alpha_size = size(BF.G_{cnt_G}{ind_alpha, ind_beta}, 1);
                    U_alpha_sigma_cell{ch_tau + 1, ch_beta + 1} = U((alpha_offset + 1) : (alpha_offset + alpha_size), :);

                    alpha_offset = alpha_offset + alpha_size;
                end
                sigma_offset = sigma_offset + sigma_size;
            end

            for ch_tau = ch_list
                ind_alpha = num_children * tau.order_ + ch_tau + 1;
                alpha = BF.tree_{ind_level_x_ch}{ind_alpha};

                G_alpha_beta = cell2mat(U_alpha_sigma_cell(ch_tau + 1, :));
                BF.G_{cnt_G}{ind_alpha, ind_beta} = G_alpha_beta;
            end
        end
    end
end

% Compress V.
level_xi = L_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});

level_x = L - level_xi;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});
R_cell = cell(m_x, m_xi);
for ind_alpha = 1 : m_x
    for ind_beta = 1 : m_xi
        [U, S, V] = MySVDSketch(BF.V_{ind_alpha, ind_beta}, tol);
        BF.V_{ind_alpha, ind_beta} = V';
        R_cell{ind_alpha, ind_beta} = U * S;
    end
end

% Compress H.
cnt_H = L_xi - h_xi + 1;
for level = (L_xi - 1) : -1 : h_xi
    level_xi = level;
    ind_level_xi = level_xi + 1;
    level_xi_ch = level_xi + 1;
    ind_level_xi_ch = level_xi_ch + 1;
    m_xi = length(BF.tree_{ind_level_xi});
    m_xi_ch = length(BF.tree_{ind_level_xi_ch});

    level_x = L - level_xi;
    ind_level_x = level_x + 1;
    level_x_par = level_x - 1;
    ind_level_x_par = level_x_par + 1;
    m_x = length(BF.tree_{ind_level_x});
    m_x_par = length(BF.tree_{ind_level_x_par});

    cnt_H = cnt_H - 1;

    % Update H.
    for ind_alpha = 1 : m_x_par
        for ind_beta = 1 : m_xi_ch
            BF.H_{cnt_H}{ind_alpha, ind_beta} = BF.H_{cnt_H}{ind_alpha, ind_beta} * R_cell{ind_alpha, ind_beta};
        end
    end

    % Compress H.
    R_cell = cell(m_x, m_xi);
    for ind_alpha = 1 : m_x_par
        for ind_sigma = 1 : m_xi
            alpha = BF.tree_{ind_level_x_par}{ind_alpha};

            sigma = BF.tree_{ind_level_xi}{ind_sigma};

            V_tau_beta_cell = cell(num_children, num_children);
            tau_offset = 0;
            for ch_alpha = ch_list
                ind_tau = num_children * alpha.order_ + ch_alpha + 1;
                tau = BF.tree_{ind_level_x}{ind_tau};

                if cnt_H - 1 == 0
                    tau_size = size(BF.M_{ind_tau, ind_sigma}, 2);
                else
                    tau_size = size(BF.H_{cnt_H - 1}{ind_tau, ind_sigma}, 2);
                end

                H_tau_beta_cell = cell(1, num_children);
                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    beta = BF.tree_{ind_level_xi_ch}{ind_beta};

                    H_alpha_beta = BF.H_{cnt_H}{ind_alpha, ind_beta};
                    H_tau_beta = H_alpha_beta((tau_offset + 1) : (tau_offset + tau_size), :);
                    H_tau_beta_cell{ch_sigma + 1} = H_tau_beta;
                end
                H_tau_sigma = cell2mat(H_tau_beta_cell);
                [U, S, V] = MySVDSketch(H_tau_sigma, tol);
                V = V';
                R_cell{ind_tau, ind_sigma} = U * S;

                beta_offset = 0;
                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    beta = BF.tree_{ind_level_xi_ch}{ind_beta};

                    beta_size = size(BF.H_{cnt_H}{ind_alpha, ind_beta}, 2);
                    V_tau_beta_cell{ch_alpha + 1, ch_sigma + 1} = V(:, (beta_offset + 1) : (beta_offset + beta_size));

                    beta_offset = beta_offset + beta_size;
                end
                tau_offset = tau_offset + tau_size;
            end

            for ch_sigma = ch_list
                ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                beta = BF.tree_{ind_level_xi_ch}{ind_beta};

                H_alpha_beta = cell2mat(V_tau_beta_cell(:, ch_sigma + 1));
                BF.H_{cnt_H}{ind_alpha, ind_beta} = H_alpha_beta;
            end
        end
    end
end

% Update M.
level_x = h_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = h_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});
for ind_tau = 1 : m_x
    for ind_sigma = 1 : m_xi
        BF.M_{ind_tau, ind_sigma} = L_cell{ind_tau, ind_sigma} * BF.M_{ind_tau, ind_sigma} * R_cell{ind_tau, ind_sigma};
    end
end

end