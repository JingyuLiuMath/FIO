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

figure_position = [100 100 1600 900];
paper_position = [0 0 24 12];
font_size = 64;
font_size_rank = font_size;
font_size_scaling = font_size;

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
    factor = 25;
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
        xlim([2^9 2^19]);
        xticks([1e3 1e4 1e5 1e6]);
        ylim([18 26]);
        yticks([18 20 22 24 26]);
    end
end
% title(title_name, "Interpreter", "latex");
lgd = legend("Location", "southeast", "Interpreter", "latex");
set(gca, 'FontSize', font_size_rank);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', paper_position);   % 统一的物理尺寸
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, figure_name + ".png", "-dpng", "-r200");
print(gcf, figure_name + ".eps", "-depsc", "-r200");

% ========== Time scaling ==========
figure();
set(gcf, 'Position', figure_position);
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
color_aBF_ref_act = "#EDB120";

color_cHSS = "#D95319";
color_cHSS_ref = "r";
color_cHSS_ref_act = "#EDB120";

color_fHSS = "#7E2F8E";
color_fHSS_ref = "m";

color_solve = "#77AC30";
color_solve_ref = "c";

apply_bf_flag = 0;

if apply_bf_flag == 1
    figure_name = figure_name + "_bf";
end

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

if apply_bf_flag == 1
    % Construct BF.
    scaling_type = "$O(N \log N)$";
    factor = mean(t_construct_BF_list);
    plot_ref_curve(N_list, scaling_type, factor, color_cBF_ref);

    % Apply BF.
    scaling_type = "$O(N \log N)$";
    factor = mean(t_apply_BF_list);
    plot_ref_curve(N_list, scaling_type, factor, color_aBF_ref);

    if isfield(result_list(1), "n")
        scaling_type = "$O(N^{1.5})$";
        factor = mean(t_apply_BF_list);
        plot_ref_curve(N_list, scaling_type, factor, color_aBF_ref_act);
    end
end

if apply_bf_flag == 0
    % Construct BF.
    scaling_type = "$O(N \log N)$";
    factor = mean(t_construct_BF_list);
    plot_ref_curve(N_list, scaling_type, factor, color_cBF_ref);

    % Construct HSS.
    if isfield(result_list(1), "n")
        scaling_type = "$O(N^{1.5} \log N)$";
        factor = mean(t_construct_HSS_list) * 0.8;
        plot_ref_curve(N_list, scaling_type, factor, color_cHSS_ref);

        scaling_type = "$O(N^{2})$";
        factor = mean(t_construct_HSS_list);
        plot_ref_curve(N_list, scaling_type, factor, color_cHSS_ref_act);
    else
        scaling_type = "$O(N \log^{2} N)$";
        factor = mean(t_construct_HSS_list);
    end

    % Factor HSS.
    if isfield(result_list(1), "n")
        scaling_type = "$O(N^{1.5})$";
    else
        scaling_type = "$O(N)$";
    end
    factor = mean(t_factor_HSS_list);
    plot_ref_curve(N_list, scaling_type, factor, color_fHSS_ref);

    % Solve HSS.
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
        if apply_bf_flag == 1
            ylim([1e-3 1e3]);
            yticks([1e-2 1e-1 1e0 1e1 1e2 1e3]);
        else
            ylim([1e-2 1e6]);
            yticks([1e-2 1e0 1e2 1e4 1e6]);
        end
    else
        xlim([2^9 2^19]);
        xticks([1e3 1e4 1e5 1e6]);
        ylim([1e-2 1e3]);
        yticks([1e-2 1e-1 1e0 1e1 1e2 1e3 1e4]);
    end
end
% title(title_name, "Interpreter", "latex");
if isfield(result_list(1), "n")
    lgd = legend("Location", "eastoutside", "Interpreter", "latex", "NumColumns", 1);
else
    lgd = legend("Location", "eastoutside", "Interpreter", "latex", "NumColumns", 1);
end
set(gca, 'FontSize', font_size_scaling);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', paper_position);   % 统一的物理尺寸
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, figure_name + ".png", "-dpng", "-r200");
print(gcf, figure_name + ".eps", "-depsc", "-r200");

end

function plot_single_curve(N_list, t_list, marker_list, display_name_list, color_list)

num_param = size(t_list, 2);

for it_param = 1 : num_param
    loglog(N_list, t_list(:, it_param), ...
        "LineWidth", 8, ...
        "Marker", marker_list(it_param), ...
        "Markersize", 16, ...
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
        linestyle = "-.";
end

loglog(N_list, ref_line, ...
    "LineWidth", 8, ...
    "LineStyle", linestyle, ...
    "DisplayName", scaling_type, ...
    "Color", color);
hold on;

end