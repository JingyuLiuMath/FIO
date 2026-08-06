function print_table_direct(result_list, caption_name, label_name)

fprintf("\n");
fprintf("\\begin{table}[tbhp]\n");
fprintf("\\centering\n");
fprintf("\\begin{tabular}{c | cc | ccc | cc}\n");
fprintf("\\toprule\n");
fprintf("\\(\\numtot\\) ");
fprintf("& \\(t_{\\construct \\BF}\\) ");
fprintf("& \\(e_{\\construct \\BF}\\) ");
fprintf("& \\(t_{\\NS}\\) ");
fprintf("& \\(n_{\\NS}\\) ");
fprintf("& \\(r_{\\NS}\\) ");
fprintf("& \\(t_{\\direct}\\) ");
fprintf("& \\(e_{\\direct}\\) ");
fprintf("\\\\ \n");

for it_n = 1 : length(result_list)
    result = result_list{it_n};
    fprintf("\\midrule\n");
    if isfield(result, "n")
        fprintf("\\(%d^{2}\\) ", result.n);
    else
        fprintf("\\(2^{%d}\\) ", log2(result.N));
    end
    fprintf("& %.1e ", result.t_construct_BF);
    fprintf("& %.1e ", result.rel_err_BF);
    fprintf("& %.1e ", result.t_construct_NS);
    fprintf("& %d ", result.history_ns.num_iterations);
    fprintf("& %.1e ", result.history_ns.residual(end));
    fprintf("& %.1e ", result.t_apply_NS);
    fprintf("& %.1e ", result.rel_err_ns);
    fprintf("\\\\ \n");
end

fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\end{table}\n");
fprintf("\n");

end
