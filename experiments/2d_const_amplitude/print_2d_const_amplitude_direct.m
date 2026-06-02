exp_2d_const_amplitude_settings;

originalPath = path;
addpath("../");

caption_name = "Results " + ... 
    "of the direct solver " + ...
    "for 2D FIO " + ...
    "with constant amplitude.";
label_name = "2d_const_amplitude_direct";

result_list = [];
for it_n = 1 : num_n
    p = p_list(it_n);

    result_list_tmp = [];
    for it_hss_tol = 1 : num_tol_hss
        tol_hss = tol_hss_list(it_hss_tol);
        load(data_path + "2d_const_amplituide_results_indep" ...
            + "_" + string(p) ...
            + "_" + string(tol_hss) ...
            + ".mat");
        result_list_tmp = [result_list_tmp, result];
    end
    result_list = [result_list; result_list_tmp];
end

if 0
    print_table_direct(result_list, ...
        caption_name, label_name, ...
        tol_hss_display_list);
else
    print_table_direct_without_tol(result_list, ...
        caption_name, label_name);
end

path(originalPath);
