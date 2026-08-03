function Construct_Interpolation(BF, a_func, phi_func, r)

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

L = BF.L_;
h_x = BF.h_x_;
h_xi = BF.h_xi_;
L_x = BF.L_x_;
L_xi = BF.L_xi_;

num_children = BF.num_children_;
ch_list = BF.ch_list_;

% Construct M.
debug_mode = 0;
level_x = h_x;
ind_level_x = level_x + 1;
m_x = length(BF.tree_{ind_level_x});

level_xi = h_xi;
ind_level_xi = level_xi + 1;
m_xi = length(BF.tree_{ind_level_xi});
BF.M_ = cell(m_x, m_xi);
if debug_mode == 1
    fprintf("  level_x: %d, level_xi: %d\n", level_x, level_xi);
end
if level_x == L_x || debug_mode == 1
    BF.U_ = cell(m_x, m_xi);
end
if level_xi == L_xi || debug_mode == 1
    BF.V_ = cell(m_x, m_xi);
end
for ind_tau = 1 : m_x
    tau = BF.tree_{ind_level_x}{ind_tau};
    % x_tau = tau.SpacePts();
    y_tau = tau.SpaceCt();
    z_tau = tau.SpaceChebPts(r);
    for ind_sigma = 1 : m_xi
        sigma = BF.tree_{ind_level_xi}{ind_sigma};
        % xi_sigma = sigma.FreqPts();
        eta_sigma = sigma.FreqCt();
        gamma_sigma = sigma.FreqChebPts(r);

        M_tau_sigma = k_func(z_tau, gamma_sigma);
        BF.M_{ind_tau, ind_sigma} = M_tau_sigma;

        if level_x == L_x || debug_mode == 1
            x_tau = tau.SpacePts();
            P = EvalLagrange(z_tau, x_tau);
            U_shift = phi_func(x_tau, eta_sigma) - phi_func(z_tau, eta_sigma).';
            U_shift = complex(cos(2 * pi * U_shift), sin(2 * pi * U_shift));
            U_tau_sigma = U_shift .* P;
            BF.U_{ind_tau, ind_sigma} = U_tau_sigma;
        end

        if level_xi == L_xi || debug_mode == 1
            xi_sigma = sigma.FreqPts();
            Q = EvalLagrange(gamma_sigma, xi_sigma).';
            V_shift = phi_func(y_tau, xi_sigma) - phi_func(y_tau, gamma_sigma).';
            V_shift = complex(cos(2 * pi * V_shift), sin(2 * pi * V_shift));
            V_tau_sigma = V_shift .* Q;
            BF.V_{ind_tau, ind_sigma} = V_tau_sigma;
        end

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
cnt_G = L_x - h_x + 1;
BF.G_ = cell(1, L_x - h_x);
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

    if debug_mode == 1
        U_old_cell = BF.U_;
        fprintf("  level_x: %d, level_xi: %d\n", level_x, level_xi);
    end

    if level_x == L_x || debug_mode == 1
        BF.U_ = cell(m_x, m_xi);
    end

    BF.G_{cnt_G} = cell(m_x, m_xi);
    for ind_alpha = 1 : m_x_par
        alpha = BF.tree_{ind_level_x_par}{ind_alpha};
        % x_alpha = alpha.SpacePts();
        % y_alpha = alpha.SpaceCt();
        z_alpha = alpha.SpaceChebPts(r);
        for ind_sigma = 1 : m_xi
            sigma = BF.tree_{ind_level_xi}{ind_sigma};
            % xi_sigma = sigma.FreqPts();
            eta_sigma = sigma.FreqCt();
            % gamma_sigma = sigma.FreqChebPts(r);

            tau_offset = 0;
            for ch_alpha = ch_list
                ind_tau = num_children * alpha.order_ + ch_alpha + 1;
                tau = BF.tree_{ind_level_x}{ind_tau};
                % x_tau = tau.SpacePts();
                % y_tau = tau.SpaceCt();
                z_tau = tau.SpaceChebPts(r);

                if level_x == L_x || debug_mode == 1
                    x_tau = tau.SpacePts();

                    P = EvalLagrange(z_tau, x_tau);
                    U_shift = phi_func(x_tau, eta_sigma) - phi_func(z_tau, eta_sigma).';
                    U_shift = complex(cos(2 * pi * U_shift), sin(2 * pi * U_shift));
                    U_tau_sigma = U_shift .* P;
                    BF.U_{ind_tau, ind_sigma} = U_tau_sigma;
                end

                P = EvalLagrange(z_alpha, z_tau);
                G_tau_beta_cell = cell(1, num_children);

                if debug_mode == 1
                    U_tau_beta_old_cell = cell(1, num_children);
                end

                for ch_sigma = ch_list
                    ind_beta = num_children * sigma.order_ + ch_sigma + 1;
                    beta = BF.tree_{ind_level_xi_ch}{ind_beta};
                    % xi_beta = beta.FreqPts();
                    eta_beta = beta.FreqCt();
                    % gamma_beta = beta.FreqChebPts(r);

                    U_shift = phi_func(z_tau, eta_beta) - phi_func(z_alpha, eta_beta).';
                    U_shift = complex(cos(2 * pi * U_shift), sin(2 * pi * U_shift));
                    G_tau_beta = U_shift .* P;
                    G_tau_beta_cell{ch_sigma + 1} = G_tau_beta;

                    if debug_mode == 1
                        U_old = U_old_cell{ind_alpha, ind_beta};
                        U_tau_beta_old_cell{ch_sigma + 1} = U_old(...
                            (tau_offset + 1) : (tau_offset + tau.ind_size_), :);
                    end
                end

                G_tau_sigma = cell2mat(G_tau_beta_cell);
                BF.G_{cnt_G}{ind_tau, ind_sigma} = G_tau_sigma;

                if debug_mode == 1
                    U_tau_sigma_old = cell2mat(U_tau_beta_old_cell);
                    rel_err_tau_sigma = norm(U_tau_sigma_old ...
                        - U_tau_sigma * G_tau_sigma, "fro") / norm(U_tau_sigma_old, "fro");
                    fprintf("  ind_tau: %d, ind_sigma: %d, rel_err: %.1e\n", ...
                        ind_tau, ind_sigma, rel_err_tau_sigma);
                end

                tau_offset = tau_offset + tau.ind_size_;
            end
        end
    end
end

% Construct H and V.
debug_mode = 0;
cnt_H = 0;
BF.H_ = cell(1, L_xi - h_xi);
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

    if debug_mode == 1
        V_old_cell = BF.V_;
        fprintf("  level_x: %d, level_xi: %d\n", level_x, level_xi);
    end

    if level_xi == L_xi || debug_mode == 1
        BF.V_ = cell(m_x, m_xi);
    end

    BF.H_{cnt_H} = cell(m_x, m_xi);
    for ind_beta = 1 : m_xi_par
        beta = BF.tree_{ind_level_xi_par}{ind_beta};
        % xi_beta = beta.FreqPts();
        % eta_beta = beta.FreqCt();
        gamma_beta = beta.FreqChebPts(r);
        for ind_tau = 1 : m_x
            tau = BF.tree_{ind_level_x}{ind_tau};
            % x_tau = tau.SpacePts();
            y_tau = tau.SpaceCt();
            % z_tau = tau.SpaceChebPts(r);

            sigma_offset = 0;
            for ch_beta = ch_list
                ind_sigma = num_children * beta.order_ + ch_beta + 1;
                sigma = BF.tree_{ind_level_xi}{ind_sigma};
                % xi_sigma = sigma.FreqPts();
                % eta_sigma = sigma.FreqCt();
                gamma_sigma = sigma.FreqChebPts(r);

                if level_xi == L_xi || debug_mode == 1
                    xi_sigma = sigma.FreqPts();

                    Q = EvalLagrange(gamma_sigma, xi_sigma).';
                    V_shift = phi_func(y_tau, xi_sigma) - phi_func(y_tau, gamma_sigma).';
                    V_shift = complex(cos(2 * pi * V_shift), sin(2 * pi * V_shift));
                    V_tau_sigma = V_shift .* Q;
                    BF.V_{ind_tau, ind_sigma} = V_tau_sigma;
                end

                Q = EvalLagrange(gamma_beta, gamma_sigma).';
                H_alpha_sigma_cell = cell(num_children, 1);

                if debug_mode == 1
                    V_alpha_sigma_old_cell = cell(num_children, 1);
                end

                for ch_tau = ch_list
                    ind_alpha = num_children * tau.order_ + ch_tau + 1;
                    alpha = BF.tree_{ind_level_x_ch}{ind_alpha};
                    % x_alpha = alpha.SpacePts();
                    y_alpha = alpha.SpaceCt();
                    % z_alpha = alpha.SpaceChebPts(r);

                    V_shift = phi_func(y_alpha, gamma_sigma) - phi_func(y_alpha, gamma_beta).';
                    V_shift = complex(cos(2 * pi * V_shift), sin(2 * pi * V_shift));
                    H_alpha_tau = V_shift .* Q;
                    H_alpha_sigma_cell{ch_tau + 1} = H_alpha_tau;

                    if debug_mode == 1
                        V_old = V_old_cell{ind_alpha, ind_beta};
                        V_alpha_sigma_old_cell{ch_tau + 1} = V_old(...
                            :, (sigma_offset + 1) : (sigma_offset + sigma.ind_size_));
                    end
                end

                H_tau_sigma = cell2mat(H_alpha_sigma_cell);
                BF.H_{cnt_H}{ind_tau, ind_sigma} = H_tau_sigma;

                if debug_mode == 1
                    V_tau_sigma_old = cell2mat(V_alpha_sigma_old_cell);
                    rel_err_tau_sigma = norm(V_tau_sigma_old ...
                        - H_tau_sigma * V_tau_sigma, "fro") / norm(V_tau_sigma_old, "fro");
                    fprintf("  ind_tau: %d, ind_sigma: %d, rel_err: %.1e\n", ...
                        ind_tau, ind_sigma, rel_err_tau_sigma);
                end

                sigma_offset = sigma_offset + sigma.ind_size_;
            end
        end
    end
end

end
