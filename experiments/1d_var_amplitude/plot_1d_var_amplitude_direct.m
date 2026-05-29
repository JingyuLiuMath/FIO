exp_1d_var_amplitude_settings;

originalPath = path;
addpath("../");

result_list = [];
for it_n = 1 : num_n
    p = p_list(it_n);

    result_list_tmp = [];
    for it_hss_tol = 1 : num_tol_hss
        tol_hss = tol_hss_list(it_hss_tol);
        load(data_path + "1d_var_amplituide_results_indep" ...
            + "_" + string(p) ...
            + "_" + string(tol_hss) ...
            + ".mat");
        result_list_tmp = [result_list_tmp, result];
    end
    result_list = [result_list; result_list_tmp];
end

plot_figure_direct(result_list, figure_prefix);

path(originalPath);
