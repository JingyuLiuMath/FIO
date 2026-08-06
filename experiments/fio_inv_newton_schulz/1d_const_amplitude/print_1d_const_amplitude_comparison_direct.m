exp_1d_const_amplitude_settings;

originalPath = path;
addpath("../");

result_bhp_list = cell(0, 1);
result_smp_list = cell(0, 1);
for it_p = 1 : num_n
    p = p_list(it_p);
    file_bhp = data_path ...
        + "1d_const_amplitude_results_newton_schulz" ...
        + "_" + string(p) + ".mat";
    file_smp = data_path ...
        + "1d_const_amplitude_results_newton_schulz_smp" ...
        + "_" + string(p) + ".mat";
    if ~isfile(file_bhp) || ~isfile(file_smp)
        continue;
    end
    result_bhp_data = load(file_bhp, "result");
    result_smp_data = load(file_smp, "result");
    result_bhp_list{end + 1, 1} = result_bhp_data.result;
    result_smp_list{end + 1, 1} = result_smp_data.result;
end

print_table_direct_comparison(...
    result_bhp_list, result_smp_list, ...
    "Comparison of the BHP- and SMP-based direct solvers " ...
    + "for a 1D FIO with constant amplitude.", ...
    "1d_const_amplitude_newton_schulz_direct_comparison");

path(originalPath);
