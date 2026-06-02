function plot_figure_direct(result_list, ...
    figure_prefix)

num_n = size(result_list, 1);

N_list = zeros(num_n, 1);
hss_rank_list = zeros(num_n, 1);

t_construct_BF_list = zeros(num_n, 1);
t_apply_BF_list = zeros(num_n, 1);
rel_err_BF_list = zeros(num_n, 1);

t_construct_HSS_list = zeros(num_n, 1);
t_factor_HSS_list = zeros(num_n, 1);
rel_err_HSS_list = zeros(num_n, 1);

t_solve_list = zeros(num_n, 1);
rel_err_direct_list = zeros(num_n, 1);

for it_n = 1 : num_n
    curr_result = result_list(it_n);

    hss_rank_list(it_n) = curr_result.hss_rank;

    t_construct_BF_list(it_n) = curr_result.t_construct_BF;
    t_apply_BF_list(it_n) = curr_result.t_apply_BF;
    rel_err_BF_list(it_n) = curr_result.rel_err_BF;

    t_construct_HSS_list(it_n) = curr_result.t_construct_HSS;
    t_factor_HSS_list(it_n) = curr_result.t_factor_HSS;
    rel_err_HSS_list(it_n) = curr_result.rel_err_HSS;

    t_solve_list(it_n) = curr_result.t_direct + curr_result.t_apply_BF;
    rel_err_direct_list(it_n) = curr_result.rel_err_direct;
    N_list(it_n) = curr_result.N;
end

% ========== Rank ==========
figure();
xlabel_name = "$N$";
ylabel_name = "numerical rank";

title_name = "Rank";
my_name = "_hss_rank";
figure_name = figure_prefix + my_name;

color_rank = "#0072BD";
color_rank_ref = "#D95319";

t_list = hss_rank_list;
marker_list = ["x"];
display_name_list = ["hss rank"];
plot_single_curve(N_list, t_list, marker_list, display_name_list, color_rank);

if isfield(result_list(1), "n")
    scaling_type = "$O(\sqrt{N})$";
    factor = mean(hss_rank_list);
else
    scaling_type = "$O(1)$";
    factor = ceil(mean(hss_rank_list) * 1.1);
end
plot_ref_curve(N_list, scaling_type, factor, color_rank_ref);

xlabel(xlabel_name, "Interpreter", "latex");
ylabel(ylabel_name, "Interpreter", "latex");
if isunix
    if isfield(result_list(1), "n")
        xlim([1e3 1e6]);
        xticks([1e3 1e4 1e5 1e6]);
        ylim([128 4096]);
        yticks([256 512 1024 2048]);
    else
        ylim([min(hss_rank_list) ceil(max(hss_rank_list) * 1.2)]);
    end
end
% title(title_name, "Interpreter", "latex");
legend("Location", "southeast", "Interpreter", "latex", "NumColumns", 2);
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");

% ========== Time scaling ==========
figure();
xlabel_name = "$N$";
ylabel_name = "time (s)";

title_name = "Time scaling";
my_name = "_time_scaling";
figure_name = figure_prefix + my_name;

marker_list = ["o", "+", "*", "x"];

color_cBF = "#0072BD";
color_cBF_ref = "b";

color_aBF = "#D95319";
color_aBF_ref = "r";

color_cHSS = "#D95319";
color_cHSS_ref = "r";

color_fHSS = "#7E2F8E";
color_fHSS_ref = "m";

color_solve = "#77AC30";
color_solve_ref = "c";

apply_bf_flag = 0;

if apply_bf_flag == 1
    t_list = [t_construct_BF_list, ...
        t_apply_BF_list];
    display_name_list = ["$t_{\mathrm{cBF}}$", ...
        "$t_{\mathrm{aBF}}$"];
    color_list = [color_cBF, ...
        color_aBF];
else
    t_list = [t_construct_BF_list, ...
        t_construct_HSS_list, ...
        t_factor_HSS_list, ...
        t_solve_list];
    display_name_list = ["$t_{\mathrm{cBF}}$", ...
        "$t_{\mathrm{cHSS}}$", ...
        "$t_{\mathrm{fHSS}}$", ...
        "$t_{\mathrm{sHSS}}$"];
    color_list = [color_cBF, ...
        color_cHSS, ...
        color_fHSS, ...
        color_solve];
end

plot_single_curve(N_list, t_list, marker_list, display_name_list, color_list);

% Construct BF.
scaling_type = "$O(N \log N)$";
factor = mean(t_construct_BF_list);
plot_ref_curve(N_list, scaling_type, factor, color_cBF_ref);

% Apply BF.
if apply_bf_flag == 1
    scaling_type = "$O(N \log N)$";
    factor = mean(t_apply_BF_list);
    plot_ref_curve(N_list, scaling_type, factor, color_aBF_ref);

    scaling_type = "$O(N^{1.5} \log N)$";
    factor = mean(t_apply_BF_list);
    plot_ref_curve(N_list, scaling_type, factor, color_aBF_ref);
end

if apply_bf_flag == 0
    % Construct HSS.
    if isfield(result_list(1), "n")
        scaling_type = "$O(N^{1.5} \log N)$";
        % scaling_type = "$O(N^{2})$";
    else
        scaling_type = "$O(N \log^{2} N)$";
    end
    factor = mean(t_construct_HSS_list);
    plot_ref_curve(N_list, scaling_type, factor, color_cHSS_ref);

    % Factor HSS.
    if isfield(result_list(1), "n")
        scaling_type = "$O(N^{1.5})$";
    else
        scaling_type = "$O(N)$";
    end
    factor = mean(t_factor_HSS_list);
    plot_ref_curve(N_list, scaling_type, factor, color_fHSS_ref);

    % Solve
    if isfield(result_list(1), "n")
        scaling_type = "$O(N \log N)$";
    else
        scaling_type = "$O(N)$";
    end
    factor = mean(t_solve_list);
    plot_ref_curve(N_list, scaling_type, factor, color_solve_ref);
end

xlabel(xlabel_name, "Interpreter", "latex");
ylabel(ylabel_name, "Interpreter", "latex");
if isunix
    if isfield(result_list(1), "n")
        xlim([1e3 1e6]);
        xticks([1e3 1e4 1e5 1e6]);
    else
        xlim([2^9 2^20])
        xticks([1e3 1e4 1e5 1e6])
        ylim([1e-4 1e3]);
        yticks([1e-3 1e-1 1e1 1e3]);
    end
end
% title(title_name, "Interpreter", "latex");

lgd = legend("Location", "southeast", "Interpreter", "latex", "NumColumns", 2);
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");

% ========== Relative error ==========
figure();
set(gca, 'Position', [0.13 0.11 0.775 0.815]);
xlabel_name = "$N$";

title_name = "Relative error";
my_name = "_rel_err";
figure_name = figure_prefix + my_name;

t_list = [rel_err_BF_list, ...
    rel_err_HSS_list, ...
    rel_err_direct_list];
marker_list = ["o", "+", "*", "x"];
display_name_list = ["$e_{\mathrm{BF}}$", ...
    "$e_{\mathrm{HSS}}$", ...
    "$e_{\mathrm{s}}$"];
color_list = ["#0072BD", "#7E2F8E", "#77AC30"];
plot_single_curve(N_list, t_list, marker_list, display_name_list, color_list);

xlabel(xlabel_name, "Interpreter", "latex");
title(title_name, "Interpreter", "latex");
legend("Location", "southeast", "Interpreter", "latex");
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");

end

function plot_single_curve(N_list, t_list, marker_list, display_name_list, color_list)

num_param = size(t_list, 2);

for it_param = 1 : num_param
    loglog(N_list, t_list(:, it_param), ...
        "LineWidth", 2, ...
        "Marker", marker_list(it_param), ...
        "Markersize", 12, ...
        "Color", color_list(it_param), ...
        "DisplayName", display_name_list(it_param));
    hold on;
end

end

function plot_ref_curve(N_list, scaling_type, factor, color)

switch scaling_type
    case ""
        return;
    case "$O(1)$"
        ref_line = ones(size(N_list));
        ref_line = ref_line / mean(ref_line) * factor;
        linestyle = "--";
    case "$O(N)$"
        ref_line = N_list;
        ref_line = ref_line / mean(ref_line) * factor;
        linestyle = ":";
    case "$O(\sqrt{N})$"
        ref_line = sqrt(N_list);
        ref_line = ref_line / mean(ref_line) * factor;
        linestyle = "--";
    case "$O(N \log N)$"
        ref_line = N_list .* log2(N_list);
        ref_line = ref_line / mean(ref_line) * factor;
        linestyle = "--";
    case "$O(N \log^{2} N)$"
        ref_line = N_list .* (log2(N_list).^2);
        ref_line = ref_line / mean(ref_line) * factor;
        linestyle = "-.";
    case "$O(N^{1.5})$"
        ref_line = N_list.^(1.5);
        ref_line = ref_line / mean(ref_line) * factor;
        linestyle = ":";
    case "$O(N^{1.5} \log N)$"
        ref_line = N_list.^(1.5) .* log2(N_list);
        ref_line = ref_line / mean(ref_line) * factor;
        linestyle = "-.";
    case "$O(N^{2})$"
        ref_line = N_list.^(2);
        ref_line = ref_line / mean(ref_line) * factor;
end

loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", linestyle, ...
    "DisplayName", scaling_type, ...
    "Color", color);
hold on;

end