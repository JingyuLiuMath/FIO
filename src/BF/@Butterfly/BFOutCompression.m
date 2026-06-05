function BFOutCompression(BF, tol)

arguments (Input)
    BF Butterfly;
    tol (1, 1) double;
end

L = BF.L_;
h_x = BF.h_x_;
h_xi = BF.h_xi_;

% Compress M.
level_x = h_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = L - level_x;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});

Lcell = cell(size(BF.M_));
Rcell = cell(size(BF.M_));
for ind_tau = 1 : m_x
    for ind_sigma = 1 : m_xi
        [U, S, V] = MySVDSketch(BF.M_{ind_tau, ind_sigma}, tol);
        Lcell{ind_tau, ind_sigma} = U*sqrt(S);
        BF.M_{ind_tau, ind_sigma} = eye(size(S));
        Rcell{ind_tau, ind_sigma} = sqrt(S)*V';
    end
end

% Compress G.
cnt_G = 0;
for level = h_x : (L - 1)
    level_x = level + 1;
    ind_level_x = level_x + 1;
    level_x_par = level_x - 1;
    ind_level_x_par = level_x_par + 1;
    m_x = length(BF.tree_{ind_level_x});

    level_xi = L - level_x;
    ind_level_xi = level_xi + 1;
    level_xi_ch = level_xi + 1;
    ind_level_xi_ch = level_xi_ch + 1;
    m_xi = length(BF.tree_{ind_level_xi});

    if cnt_G == length(BF.G_)
        level_x = level_x - 1;
        ind_level_x = level_x + 1;
        break;
    end
    cnt_G = cnt_G + 1;
    % Update G.
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            tau = BF.tree_{ind_level_x}{ind_tau};

            sigma = BF.tree_{ind_level_xi}{ind_sigma};

            alpha_order = floor(tau.order_ / 2);
            ind_alpha = alpha_order + 1;
            alpha = BF.tree_{ind_level_x_par}{ind_alpha};

            L_tau_sigma_cell = cell(1, 2);
            for cind_beta = [0, 1]
                ind_beta = 2 * sigma.order_ + cind_beta + 1;
                beta = BF.tree_{ind_level_xi_ch}{ind_beta};

                L_alpha_beta = Lcell{ind_alpha, ind_beta};
                offset = 0;
                for cind_alpha = [0, 1]
                    order_ch_alpha = 2 * alpha.order_ + cind_alpha;
                    ind_ch_alpha = order_ch_alpha + 1;
                    if order_ch_alpha == tau.order_
                        break;
                    end
                    offset = offset + BF.tree_{ind_level_x}{ind_ch_alpha}.ind_size_;
                end
                L_tau_sigma_cell{cind_beta + 1} = L_alpha_beta((offset + 1) : (offset + tau.ind_size_), :);
            end

            L_tau_sigma = cell2mat(L_tau_sigma_cell);
            BF.G_{cnt_G}{ind_tau, ind_sigma} = BF.G_{cnt_G}{ind_tau, ind_sigma} * L_tau_sigma;
        end
    end

    % Compress G.
    Lcell = cell(m_x, m_xi);
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            [U, S, V] = MySVDSketch(BF.G_{cnt_G}{ind_tau, ind_sigma}, tol);
            Lcell{ind_tau, ind_sigma} = U * S;
            BF.G_{cnt_G}{ind_tau, ind_sigma} = V';
        end
    end
end

% Update U.
for ind_tau = 1 : size(BF.U_, 1)
    for ind_sigma = 1 : size(BF.U_, 2)
        BF.U_{ind_tau, ind_sigma} = BF.U_{ind_tau, ind_sigma} * Lcell{ind_tau, ind_sigma};
    end
end

% Compress H.
cnt_H = 0;
for level = h_xi : (L - 1)
    level_xi = level + 1;
    ind_level_xi = level_xi + 1;
    level_xi_par = level_xi - 1;
    ind_level_xi_par = level_xi_par + 1;
    m_xi = length(BF.tree_{ind_level_xi});

    level_x = L - level_xi;
    ind_level_x = level_x + 1;
    level_x_ch = level_x + 1;
    ind_level_x_ch = level_x_ch + 1;
    m_x = length(BF.tree_{ind_level_x});

    if cnt_H == length(BF.H_)
        level_xi = level_xi - 1;
        ind_level_xi = level_xi + 1;
        break;
    end
    cnt_H = cnt_H + 1;

    % Update H.
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            tau = BF.tree_{ind_level_x}{ind_tau};

            sigma = BF.tree_{ind_level_xi}{ind_sigma};

            beta_order = floor(sigma.order_ / 2);
            ind_beta = beta_order + 1;
            beta = BF.tree_{ind_level_xi_par}{ind_beta};

            R_tau_sigma_cell = cell(2, 1);
            for cind_alpha = [0, 1]
                ind_alpha = 2 * tau.order_ + cind_alpha + 1;
                alpha = BF.tree_{ind_level_x_ch}{ind_alpha};

                R_alpha_beta = Rcell{ind_alpha, ind_beta};
                offset = 0;
                for cind_beta = [0, 1]
                    order_ch_beta = 2 * beta.order_ + cind_beta;
                    ind_ch_beta = order_ch_beta + 1;
                    if order_ch_beta == sigma.order_
                        break;
                    end
                    offset = offset + BF.tree_{ind_level_xi}{ind_ch_beta}.ind_size_;
                end
                R_tau_sigma_cell{cind_alpha + 1} = R_alpha_beta(:, (offset + 1) : (offset + sigma.ind_size_));
            end

            R_tau_sigma = cell2mat(R_tau_sigma_cell);
            BF.H_{cnt_H}{ind_tau, ind_sigma} = R_tau_sigma * BF.H_{cnt_H}{ind_tau, ind_sigma};
        end
    end

    % Compress H.
    Rcell = cell(m_x, m_xi);
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            [U, S, V] = MySVDSketch(BF.H_{cnt_H}{ind_tau, ind_sigma}, tol);
            Lcell{ind_tau, ind_sigma} = S * V';
            BF.G_{cnt_G}{ind_tau, ind_sigma} = U;
        end
    end
end

% Update V.
for ind_tau = 1 : size(BF.V_, 1)
    for ind_sigma = 1 : size(BF.V_, 2)
        BF.V_{ind_tau, ind_sigma} = Rcell{ind_tau, ind_sigma} * BF.V_{ind_tau, ind_sigma};
    end
end

end