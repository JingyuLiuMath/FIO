function print_table_direct_comparison(...
    result_bhp_list, result_smp_list, caption_name, label_name)

assert(length(result_bhp_list) == length(result_smp_list));

fprintf("\n");
fprintf("\\begin{table}[tbhp]\n");
fprintf("\\centering\n");
fprintf("\\begin{tabular}{cc | cc | ccc | cc}\n");
fprintf("\\toprule\n");
fprintf("\\(\\numtot\\) ");
fprintf("& Form ");
fprintf("& \\(t_{\\construct \\BF}\\) ");
fprintf("& \\(e_{\\construct \\BF}\\) ");
fprintf("& \\(t_{\\NS}\\) ");
fprintf("& \\(n_{\\NS}\\) ");
fprintf("& \\(r_{\\NS}\\) ");
fprintf("& \\(t_{\\direct}\\) ");
fprintf("& \\(e_{\\direct}\\) ");
fprintf("\\\\ \n");

for it_n = 1 : length(result_bhp_list)
    result_bhp = result_bhp_list{it_n};
    result_smp = result_smp_list{it_n};
    assert(result_bhp.N == result_smp.N);

    fprintf("\\midrule\n");
    PrintRow(result_bhp, "BHP");
    PrintRow(result_smp, "SMP");
end

fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\end{table}\n");
fprintf("\n");

end

function PrintRow(result, form)

fprintf("\\(2^{%d}\\) ", log2(result.N));
fprintf("& %s ", form);
fprintf("& %.1e ", result.t_construct_BF);
fprintf("& %.1e ", result.rel_err_BF);
fprintf("& %.1e ", result.t_construct_NS);
fprintf("& %d ", result.history_ns.num_iterations);
fprintf("& %.1e ", result.history_ns.residual(end));
fprintf("& %.1e ", result.t_apply_NS);
fprintf("& %.1e ", result.rel_err_ns);
fprintf("\\\\ \n");

end
