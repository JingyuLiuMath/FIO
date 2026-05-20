function plot_figure_direct(result_list, ...
    figure_prefix)

num_n = size(result_list, 1);

N_list = zeros(num_n, 1);
hss_rank_list = zeros(num_n, 1);

t_construct_BF_list = zeros(num_n, 1);
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
    rel_err_BF_list(it_n) = curr_result.rel_err_BF;

    t_construct_HSS_list(it_n) = curr_result.t_construct_HSS;
    t_factor_HSS_list(it_n) = curr_result.t_factor_HSS;
    rel_err_HSS_list(it_n) = curr_result.rel_err_HSS;

    t_solve_list(it_n) = curr_result.t_direct + curr_result.t_apply_BF;
    rel_err_direct_list(it_n) = curr_result.rel_err_direct;
    N_list(it_n) = curr_result.N;
end

figure();
xlabel_name = "$N$";

title_name = "Rank";
my_name = "_hss_rank";
figure_name = figure_prefix + my_name;

t_list = hss_rank_list;
marker_list = ["x"];
display_name_list = ["hss rank"];
plot_single_curve(N_list, t_list, marker_list, display_name_list);

if isfield(result_list(1), "n")
    scaling_type = "$O(\sqrt{N})$";
else
    scaling_type = "$O(1)$";
end
factor = hss_rank_list(end) * 1.5;
plot_ref_curve(N_list, scaling_type, factor);

xlabel(xlabel_name, "Interpreter", "latex");
title(title_name, "Interpreter", "latex");
legend("Location", "southeast", "Interpreter", "latex");
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");

figure();
xlabel_name = "$N$";
ylabel_name = "time (s)";

title_name = "Time scaling";
my_name = "_time_scaling";
figure_name = figure_prefix + my_name;

t_list = [t_construct_BF_list, ...
    t_construct_HSS_list, ...
    t_factor_HSS_list, ...
    t_solve_list];
marker_list = ["o", "+", "*", "x"];
display_name_list = ["$t_{\mathrm{cBF}}$", ...
    "$t_{\mathrm{cHSS}}$", ...
    "$t_{\mathrm{fHSS}}$", ...
    "$t_{\mathrm{s}}$"];
plot_single_curve(N_list, t_list, marker_list, display_name_list);

% Construct BF.
scaling_type = "$O(N \log (N))$";
factor = t_construct_BF_list(end) * 1.5;
plot_ref_curve(N_list, scaling_type, factor);

% Construct HSS.
if isfield(result_list(1), "n")
    scaling_type = "$O(N^{1.5} \log (N)$";
else
    scaling_type = "$O(N \log^{2}(N))$";
end
factor = t_construct_HSS_list(end) * 1.5;
plot_ref_curve(N_list, scaling_type, factor);

% Factor HSS.
if isfield(result_list(1), "n")
    scaling_type = "$O(N^{1.5})$";
else
    scaling_type = "$O(N)$";
end
factor = t_factor_HSS_list(end) * 1.5;
plot_ref_curve(N_list, scaling_type, factor);

% Solve
if isfield(result_list(1), "n")
    scaling_type = "";
else
    scaling_type = "";
end
factor = t_solve_list(end) * 1.5;
plot_ref_curve(N_list, scaling_type, factor);


xlabel(xlabel_name, "Interpreter", "latex");
ylabel(ylabel_name, "Interpreter", "latex");
title(title_name, "Interpreter", "latex");
legend("Location", "southeast", "Interpreter", "latex");
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");


figure();
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
plot_single_curve(N_list, t_list, marker_list, display_name_list);

xlabel(xlabel_name, "Interpreter", "latex");
title(title_name, "Interpreter", "latex");
legend("Location", "southeast", "Interpreter", "latex");
set(gca, 'FontSize', 24);
saveas(gcf, figure_name + ".png", "png");
saveas(gcf, figure_name + ".eps", "epsc");

end

function plot_single_curve(N_list, t_list, marker_list, display_name_list)

num_param = size(t_list, 2);

for it_param = 1 : num_param
    loglog(N_list, t_list(:, it_param), ...
        "LineWidth", 2, ...
        "Marker", marker_list(it_param), ...
        "Markersize", 12, ...
        "DisplayName", display_name_list(it_param));
    hold on;
end

end

function plot_ref_curve(N_list, scaling_type, factor)

switch scaling_type
    case ""
        return;
    case "$O(1)$"
        ref_line = ones(size(N_list));
        ref_line = ref_line / ref_line(1) * factor;
    case "$O(N)$"
        ref_line = N_list;
        ref_line = ref_line / ref_line(1) * factor;
    case "$O(\sqrt{N})$"
        ref_line = sqrt(N_list);
        ref_line = ref_line / ref_line(1) * factor;
    case "$O(N \log (N))$"
        ref_line = N_list .* log2(N_list);
        ref_line = ref_line / ref_line(1) * factor;
    case "$O(N \log^{2}(N))$"
        ref_line = N_list .* (log2(N_list).^2);
        ref_line = ref_line / ref_line(1) * factor;
    case "$O(N^{1.5})$"
        ref_line = N_list.^(1.5);
        ref_line = ref_line / ref_line(1) * factor;
    case "$O(N^{1.5} \log (N)$"
        ref_line = N_list.^(1.5) .* log2(N_list);
        ref_line = ref_line / ref_line(1) * factor;
end

loglog(N_list, ref_line, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "DisplayName", scaling_type);
hold on;

end