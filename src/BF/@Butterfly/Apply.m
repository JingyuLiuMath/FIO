function f = Apply(BF, f)

arguments (Input)
    BF Butterfly;
    f (: ,:) double;
end

arguments (Output)
    f (:, :) double;
end

N = BF.N_;
L = BF.L_;
h_x = BF.h_x_;
h_xi = BF.h_xi_;
L_x = BF.L_x_;
L_xi = BF.L_xi_;

num_children = BF.num_children_;
ch_list = BF.ch_list_;

num_col = size(f, 2);

% Apply V.
level_xi = L_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});

level_x = L - level_xi;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});
f_cell = cell(m_x, m_xi);
for ind_tau = 1 : m_x
    sigma_offset = 0;
    for ind_sigma = 1 : m_xi
        sigma_size = size(BF.V_{ind_tau, ind_sigma}, 2);
        f_cell{ind_tau, ind_sigma} = BF.V_{ind_tau, ind_sigma} * f((sigma_offset + 1) : (sigma_offset + sigma_size), :);

        sigma_offset = sigma_offset + sigma_size;
    end
end

% Apply H.
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

    % Update f.
    for ind_alpha = 1 : m_x_par
        for ind_beta = 1 : m_xi_ch
            f_cell{ind_alpha, ind_beta} = BF.H_{cnt_H}{ind_alpha, ind_beta} * f_cell{ind_alpha, ind_beta};
        end
    end

    % Reshape f.
    g_cell = cell(m_x, m_xi);
    for ind_alpha = 1 : m_x_par
        for ind_sigma = 1 : m_xi
            alpha = BF.tree_{ind_level_x_par}{ind_alpha};

            sigma = BF.tree_{ind_level_xi}{ind_sigma};

            tau_offset = 0;
            for ch_alpha = ch_list
                ind_tau = num_children * alpha.order_ + ch_alpha + 1;
                tau = BF.tree_{ind_level_x}{ind_tau};

                if cnt_H - 1 >= 1
                    tau_size = size(BF.H_{cnt_H - 1}{ind_tau, ind_sigma}, 2);
                else
                    tau_size = size(BF.M_{ind_tau, ind_sigma}, 2);
                end

                g_tau_sigma = 0;
                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    beta = BF.tree_{ind_level_xi_ch}{ind_beta};

                    f_alpha_beta = f_cell{ind_alpha, ind_beta};
                    g_tau_sigma = g_tau_sigma + f_alpha_beta((tau_offset + 1) : (tau_offset + tau_size), :);
                end
                g_cell{ind_tau, ind_sigma} = g_tau_sigma;

                tau_offset = tau_offset + tau_size;
            end
        end
    end
    f_cell = g_cell;
end

% Apply M.
level_x = h_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = h_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});
for ind_tau = 1 : m_x
    for ind_sigma = 1 : m_xi
        f_cell{ind_tau, ind_sigma} = BF.M_{ind_tau, ind_sigma} * f_cell{ind_tau, ind_sigma};
    end
end

% Apply G.
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

    % Reshape f.
    g_cell = cell(m_x, m_xi);
    for ind_alpha = 1 : m_x_par
        for ind_sigma = 1 : m_xi
            alpha = BF.tree_{ind_level_x_par}{ind_alpha};

            sigma = BF.tree_{ind_level_xi}{ind_sigma};

            for ch_alpha = ch_list
                ind_tau = num_children * alpha.order_ + ch_alpha + 1;
                tau = BF.tree_{ind_level_x}{ind_tau};

                tau_size = size(BF.G_{cnt_G}{ind_tau, ind_sigma}, 2);

                g_tau_sigma = zeros(tau_size, num_col);
                beta_offset = 0;
                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    beta = BF.tree_{ind_level_xi_ch}{ind_beta};

                    if cnt_G + 1 <= length(BF.G_)
                        beta_size = size(BF.G_{cnt_G + 1}{ind_alpha, ind_beta}, 1);
                    else
                        beta_size = size(BF.M_{ind_alpha, ind_beta}, 1);
                    end
                    f_alpha_beta = f_cell{ind_alpha, ind_beta};
                    g_tau_sigma((beta_offset + 1) : (beta_offset + beta_size), :) = f_alpha_beta;
                    
                    beta_offset = beta_offset + beta_size;
                end
                g_cell{ind_tau, ind_sigma} = g_tau_sigma;
            end
        end
    end
    f_cell = g_cell;

    % Update f.
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            f_cell{ind_tau, ind_sigma} = BF.G_{cnt_G}{ind_tau, ind_sigma} * f_cell{ind_tau, ind_sigma};
        end
    end
end

% Apply U.
level_x = L_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = L - level_x;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});
f = zeros(N, num_col);
for ind_sigma = 1 : m_xi
    tau_offset = 0;
    for ind_tau = 1 : m_x
        tau_size = size(BF.U_{ind_tau, ind_sigma}, 1);
        f((tau_offset + 1) : (tau_offset+ tau_size), :) = f((tau_offset + 1) : (tau_offset+ tau_size), :) + BF.U_{ind_tau, ind_sigma} * f_cell{ind_tau, ind_sigma};

        tau_offset = tau_offset + tau_size;
    end
end

end