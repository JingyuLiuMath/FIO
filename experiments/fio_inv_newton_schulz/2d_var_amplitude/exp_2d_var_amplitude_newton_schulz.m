exp_2d_var_amplitude_settings;

for it_p = 1 : num_n
    p = p_list(it_p);

    fprintf("\n\n\n\n");
    fprintf("p: %d\n", p);

    result_data = load(data_path ...
        + "2d_var_amplitude_results_gmres" ...
        + "_" + string(p) + ".mat", "result");
    result = run_FIO2D_inv_NewtonSchulz( ...
        result_data.result, r_ns, tol_ns, maxit_ns, num_test, ...
        type_bf_inv);

    save(data_path + "2d_var_amplitude_results_newton_schulz" ...
        + "_" + string(p) + ".mat", "result", "-v7.3");
end
