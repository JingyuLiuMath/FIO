function f = ApplyAdjoint(BF, f)

arguments (Input)
    BF Butterfly2D;
    f (: ,:) double;
end

arguments (Output)
    f (:, :) double;
end

assert(BF.constructed_, "Butterfly2D:ApplyAdjoint:NotConstructed", ...
    "Call Construct or BlackBoxConstruct before ApplyAdjoint.");

N = BF.N_;
L = BF.L_;
h_x = BF.h_x_;
h_xi = BF.h_xi_;
L_x = BF.L_x_;
L_xi = BF.L_xi_;

num_children = BF.num_children_;
ch_list = BF.ch_list_;

num_col = size(f, 2);

% Perm.
f = f(BF.x_perm_, :);

% Apply U.
level_x = L_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = L - level_x;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});
f_cell = cell(m_xi, m_x);
for ind_sigma = 1 : m_xi
    tau_offset = 0;
    for ind_tau = 1 : m_x
        tau_size = size(BF.U_{ind_tau, ind_sigma}, 1);

        f_cell{ind_sigma, ind_tau} = BF.U_{ind_tau, ind_sigma}' * f((tau_offset + 1) : (tau_offset + tau_size), :);

        tau_offset = tau_offset + tau_size;
    end
end

% Apply G.
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

    % Update f.
    for ind_beta = 1 : m_xi_par
        for ind_alpha = 1 : m_x_ch
            f_cell{ind_beta, ind_alpha} = BF.G_{cnt_G}{ind_alpha, ind_beta}' * f_cell{ind_beta, ind_alpha};
        end
    end

    % Reshape f.
    g_cell = cell(m_xi, m_x);
    for ind_beta = 1 : m_xi_par
        for ind_tau = 1 : m_x
            beta = BF.tree_{ind_level_xi_par}{ind_beta};

            tau = BF.tree_{ind_level_x}{ind_tau};

            sigma_offset = 0;
            for ch_beta = ch_list
                ind_sigma = num_children * beta.order_ + ch_beta + 1;
                sigma = BF.tree_{ind_level_xi}{ind_sigma};

                if cnt_G + 1 <= length(BF.G_)
                    sigma_size = size(BF.G_{cnt_G + 1}{ind_tau, ind_sigma}, 1);
                else
                    sigma_size = size(BF.M_{ind_tau, ind_sigma}, 1);
                end

                g_sigma_tau = 0;
                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    alpha = BF.tree_{ind_level_x_ch}{ind_alpha};

                    f_beta_alpha = f_cell{ind_beta, ind_alpha};
                    g_sigma_tau = g_sigma_tau + f_beta_alpha((sigma_offset + 1) : (sigma_offset + sigma_size), :);
                end
                g_cell{ind_sigma, ind_tau} = g_sigma_tau;

                sigma_offset = sigma_offset + sigma_size;
            end
        end
    end
    f_cell = g_cell;
end

% Apply M.
level_xi = h_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});

level_x = h_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});
for ind_sigma = 1 : m_xi
    for ind_tau = 1 : m_x
        f_cell{ind_sigma, ind_tau} = BF.M_{ind_tau, ind_sigma}' * f_cell{ind_sigma, ind_tau};
    end
end

% Apply H.
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

    % Reshape f.
    g_cell = cell(m_xi, m_x);
    for ind_beta = 1 : m_xi_par
        for ind_tau = 1 : m_x
            beta = BF.tree_{ind_level_xi_par}{ind_beta};

            tau = BF.tree_{ind_level_x}{ind_tau};

            for ch_beta = ch_list
                ind_sigma = num_children * beta.order_ + ch_beta + 1;
                sigma = BF.tree_{ind_level_xi}{ind_sigma};

                sigma_size = size(BF.H_{cnt_H}{ind_tau, ind_sigma}, 1);

                g_sigma_tau = zeros(sigma_size, num_col);
                alpha_offset = 0;
                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    alpha = BF.tree_{ind_level_x_ch}{ind_alpha};

                    if cnt_H - 1 >= 1
                        alpha_size = size(BF.H_{cnt_H - 1}{ind_alpha, ind_beta}, 2);
                    else
                        alpha_size = size(BF.M_{ind_alpha, ind_beta}, 2);
                    end
                    f_beta_alpha = f_cell{ind_beta, ind_alpha};
                    g_sigma_tau((alpha_offset + 1) : (alpha_offset + alpha_size), :) = f_beta_alpha;
                    
                    alpha_offset = alpha_offset + alpha_size;
                end
                g_cell{ind_sigma, ind_tau} = g_sigma_tau;
            end
        end
    end
    f_cell = g_cell;

    % Update f.
    for ind_sigma = 1 : m_xi
        for ind_tau = 1 : m_x
            f_cell{ind_sigma, ind_tau} = BF.H_{cnt_H}{ind_tau, ind_sigma}' * f_cell{ind_sigma, ind_tau};
        end
    end
end

% Apply V.
level_xi = L_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});

level_x = L - level_xi;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});
f = zeros(N, num_col);
for ind_tau = 1 : m_x
    sigma_offset = 0;
    for ind_sigma = 1 : m_xi
        sigma_size = size(BF.V_{ind_tau, ind_sigma}, 2);
        f((sigma_offset + 1) : (sigma_offset+ sigma_size), :) = f((sigma_offset + 1) : (sigma_offset+ sigma_size), :) + BF.V_{ind_tau, ind_sigma}' * f_cell{ind_sigma, ind_tau};

        sigma_offset = sigma_offset + sigma_size;
    end
end

% Perm.
f = f(BF.xi_perm_inv_, :);


end
