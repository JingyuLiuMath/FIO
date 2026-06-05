function Construct(BF, a_func, phi_func, r)

arguments (Input)
    BF Butterfly;
    a_func function_handle;
    phi_func function_handle;
    r (1, 1) double;
end

exp_phi_func = @(x, xi) complex(...
    cos(2 * pi * phi_func(x, xi)), ...
    sin(2 * pi * phi_func(x, xi)));
k_func = @(x, xi) a_func(x, xi) .* exp_phi_func(x, xi);

min_points = 16;

N = BF.N_;

L = BF.L_;
h_x = BF.h_x_;
h_xi = BF.h_xi_;

% Construct M.
debug_mode = 0;
level_x = h_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = h_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});

fprintf("  level_x: %d, level_xi: %d\n", level_x, level_xi);
BF.M_ = cell(m_x, m_xi);
BF.U_ = cell(m_x, m_xi);
BF.V_ = cell(m_x, m_xi);
for ind_tau = 1 : m_x
    for ind_sigma = 1 : m_xi
        tau = BF.tree_{ind_level_x}{ind_tau};
        x_tau = tau.SpacePts();
        y_tau = tau.SpaceCt();
        z_tau = tau.SpaceChebPts(r);

        sigma = BF.tree_{ind_level_xi}{ind_sigma};
        xi_sigma = sigma.FreqPts();
        eta_sigma = sigma.FreqCt();
        gamma_sigma = sigma.FreqChebPts(r);

        M_tau_sigma = k_func(z_tau, gamma_sigma);

        P = EvalLagrange(z_tau, x_tau);
        U_shift = phi_func(x_tau, eta_sigma) - phi_func(z_tau, eta_sigma).';
        U_shift = complex(cos(2 * pi * U_shift), sin(2 * pi * U_shift));
        U_tau_sigma = U_shift .* P;

        Q = EvalLagrange(gamma_sigma, xi_sigma).';
        V_shift = phi_func(y_tau, xi_sigma) - phi_func(y_tau, gamma_sigma).';
        V_shift = complex(cos(2 * pi * V_shift), sin(2 * pi * V_shift));
        V_tau_sigma = V_shift .* Q;

        BF.U_{ind_tau, ind_sigma} = U_tau_sigma;
        BF.M_{ind_tau, ind_sigma} = M_tau_sigma;
        BF.V_{ind_tau, ind_sigma} = V_tau_sigma;

        if debug_mode == 1
            K_tau_sigma = k_func(x_tau, xi_sigma);
            rel_err_tau_sigma = norm(K_tau_sigma ...
                - U_tau_sigma * M_tau_sigma * V_tau_sigma, "fro") / norm(K_tau_sigma, "fro");
            fprintf("  ind_tau: %d, ind_sigma: %d, rel_err: %.1e\n", ...
                ind_tau, ind_sigma, rel_err_tau_sigma);
        end
    end
end

% Construct G and U.
debug_mode = 0;
Gcell = cell(1, L - h_x);
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

    if N / 2^level_x <= min_points
        level_x = level_x - 1;
        ind_level_x = level_x + 1;
        break;
    end
    cnt_G = cnt_G + 1; 

    fprintf("  level_x: %d, level_xi: %d\n", level_x, level_xi);
    if debug_mode == 1
        U_old_cell = BF.U_;
    end
    BF.U_ = cell(m_x, m_xi);
    Gcell{cnt_G} = cell(m_x, m_xi);
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            tau = BF.tree_{ind_level_x}{ind_tau};
            x_tau = tau.SpacePts();
            y_tau = tau.SpaceCt();
            z_tau = tau.SpaceChebPts(r);

            sigma = BF.tree_{ind_level_xi}{ind_sigma};
            xi_sigma = sigma.FreqPts();
            eta_sigma = sigma.FreqCt();
            gamma_sigma = sigma.FreqChebPts(r);

            alpha_order = floor(tau.order_ / 2);
            ind_alpha = alpha_order + 1;
            alpha = BF.tree_{ind_level_x_par}{ind_alpha};
            x_alpha = alpha.SpacePts();
            y_alpha = alpha.SpaceCt();
            z_alpha = alpha.SpaceChebPts(r);

            P = EvalLagrange(z_tau, x_tau);
            U_shift = phi_func(x_tau, eta_sigma) - phi_func(z_tau, eta_sigma).';
            U_shift = complex(cos(2 * pi * U_shift), sin(2 * pi * U_shift));
            U_tau_sigma = U_shift .* P;

            if debug_mode == 1
                U_tau_sigma_old_cell = cell(1, 2);
            end
            G_tau_sigma_cell = cell(1, 2);
            for cind_beta = [0, 1]
                ind_beta = 2 * sigma.order_ + cind_beta + 1;
                beta = BF.tree_{ind_level_xi_ch}{ind_beta};
                xi_beta = beta.FreqPts();
                eta_beta = beta.FreqCt();
                gamma_beta = beta.FreqChebPts(r);

                P = EvalLagrange(z_alpha, z_tau);
                U_shift = phi_func(z_tau, eta_beta) - phi_func(z_alpha, eta_beta).';
                U_shift = complex(cos(2 * pi * U_shift), sin(2 * pi * U_shift));
                G_tau_beta = U_shift .* P;

                G_tau_sigma_cell{cind_beta + 1} = G_tau_beta;
                if debug_mode == 1
                    U_old = U_old_cell{ind_alpha, ind_beta};
                    offset = 0;
                    for cind_alpha = [0, 1]
                        order_ch_alpha = 2 * alpha.order_ + cind_alpha;
                        ind_ch_alpha = order_ch_alpha + 1;
                        if order_ch_alpha == tau.order_
                            break;
                        end
                        offset = offset + BF.tree_{ind_level_x}{ind_ch_alpha}.ind_size_;
                    end
                    U_tau_sigma_old_cell{cind_beta + 1} = U_old((offset + 1) : (offset + tau.ind_size_), :);
                end
            end
            G_tau_sigma = cell2mat(G_tau_sigma_cell);

            BF.U_{ind_tau, ind_sigma} = U_tau_sigma;
            Gcell{cnt_G}{ind_tau, ind_sigma} = G_tau_sigma;

            if debug_mode == 1
                U_tau_sigma_old = cell2mat(U_tau_sigma_old_cell);
                rel_err_tau_sigma = norm(U_tau_sigma_old ...
                    - U_tau_sigma * G_tau_sigma, "fro") / norm(U_tau_sigma_old, "fro");
                fprintf("  ind_tau: %d, ind_sigma: %d, rel_err: %.1e\n", ...
                    ind_tau, ind_sigma, rel_err_tau_sigma);
            end
        end
    end
end
BF.G_ = cell(1, cnt_G);
for t = 1 : cnt_G
    BF.G_{t} = Gcell{t};
end

% Construct H and V.
debug_mode = 0;
Hcell = cell(1, L - h_xi);
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

    if N / 2^level_xi <= min_points
        level_xi = level_xi - 1;
        ind_level_xi = level_xi + 1;
        break;
    end
    cnt_H = cnt_H + 1; 

    fprintf("  level_x: %d, level_xi: %d\n", level_x, level_xi);
    if debug_mode == 1
        V_old_cell = BF.V_;
    end
    BF.V_ = cell(m_x, m_xi);
    Hcell{cnt_H} = cell(m_x, m_xi);
    for ind_tau = 1 : m_x
        for ind_sigma = 1 : m_xi
            tau = BF.tree_{ind_level_x}{ind_tau};
            x_tau = tau.SpacePts();
            y_tau = tau.SpaceCt();
            z_tau = tau.SpaceChebPts(r);

            sigma = BF.tree_{ind_level_xi}{ind_sigma};
            xi_sigma = sigma.FreqPts();
            eta_sigma = sigma.FreqCt();
            gamma_sigma = sigma.FreqChebPts(r);

            beta_order = floor(sigma.order_ / 2);
            ind_beta = beta_order + 1;
            beta = BF.tree_{ind_level_xi_par}{ind_beta};
            xi_beta = beta.FreqPts();
            eta_beta = beta.FreqCt();
            gamma_beta = beta.FreqChebPts(r);

            Q = EvalLagrange(gamma_sigma, xi_sigma).';
            V_shift = phi_func(y_tau, xi_sigma) - phi_func(y_tau, gamma_sigma).';
            V_shift = complex(cos(2 * pi * V_shift), sin(2 * pi * V_shift));
            V_tau_sigma = V_shift .* Q;

            if debug_mode == 1
                V_tau_sigma_old_cell = cell(2, 1);
            end
            H_alpha_tau_cell = cell(2, 1);
            for cind_alpha = [0, 1]
                ind_alpha = 2 * tau.order_ + cind_alpha + 1;
                alpha = BF.tree_{ind_level_x_ch}{ind_alpha};
                x_alpha = alpha.SpacePts();
                y_alpha = alpha.SpaceCt();
                z_alpha = alpha.SpaceChebPts(r);

                Q = EvalLagrange(gamma_beta, gamma_sigma).';
                V_shift = phi_func(y_alpha, gamma_sigma) - phi_func(y_alpha, gamma_beta).';
                V_shift = complex(cos(2 * pi * V_shift), sin(2 * pi * V_shift));
                H_alpha_tau = V_shift .* Q;

                H_alpha_tau_cell{cind_alpha + 1} = H_alpha_tau;
                if debug_mode == 1
                    V_old = V_old_cell{ind_alpha, ind_beta};
                    offset = 0;
                    for cind_beta = [0, 1]
                        order_ch_beta = 2 * beta.order_ + cind_beta;
                        ind_ch_beta = order_ch_beta + 1;
                        if order_ch_beta == sigma.order_
                            break;
                        end
                        offset = offset + BF.tree_{ind_level_xi}{ind_ch_beta}.ind_size_;
                    end
                    V_tau_sigma_old_cell{cind_alpha + 1} = V_old(:, (offset + 1) : (offset + sigma.ind_size_));
                end
            end
            H_tau_sigma = cell2mat(H_alpha_tau_cell);

            BF.V_{ind_tau, ind_sigma} = V_tau_sigma;
            Hcell{cnt_H}{ind_tau, ind_sigma} = H_tau_sigma;

            if debug_mode == 1
                V_tau_sigma_old = cell2mat(V_tau_sigma_old_cell);
                rel_err_tau_sigma = norm(V_tau_sigma_old ...
                    - H_tau_sigma * V_tau_sigma, "fro") / norm(V_tau_sigma_old, "fro");
                fprintf("  ind_tau: %d, ind_sigma: %d, rel_err: %.1e\n", ...
                    ind_tau, ind_sigma, rel_err_tau_sigma);
            end

        end
    end
end
BF.H_ = cell(1, cnt_H);
for t = 1 : cnt_H
    BF.H_{t} = Hcell{t};
end

end