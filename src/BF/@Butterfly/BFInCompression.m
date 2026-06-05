function BFInCompression(BF, tol)

arguments (Input)
    BF Butterfly;
    tol (1, 1) double;
end

% Compress U.
Lcell = cell(size(BF.U_));
for ind_tau = 1 : size(BF.U_, 1)
    for ind_sigma = 1 : size(BF.U_, 2)
        [U, S, V] = MySVDSketch(BF.U_{ind_tau, ind_sigma}, tol);
        BF.U_{ind_tau, ind_sigma} = U;
        Lcell{ind_tau, ind_sigma} = S * V';
    end
end

% Compress G.
cnt_G = length(BF.G_);
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


end