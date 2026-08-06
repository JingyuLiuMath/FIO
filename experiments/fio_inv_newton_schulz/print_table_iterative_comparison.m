function print_table_iterative_comparison(...
    result_bhp_list, result_smp_list, caption_name, label_name)

assert(length(result_bhp_list) == length(result_smp_list));

fprintf("\n");
fprintf("\\begin{table}[tbhp]\n");
fprintf("\\centering\n");
fprintf("\\begin{tabular}{cc | ccc | cccc}\n");
fprintf("\\toprule\n");
fprintf("\\(\\numtot\\) ");
fprintf("& Form ");
fprintf("& \\(t_{\\iter}\\) ");
fprintf("& \\(n_{\\iter}\\) ");
fprintf("& \\(e_{\\iter}\\) ");
fprintf("& \\(t_{\\offline}\\) ");
fprintf("& \\(t_{\\piter}\\) ");
fprintf("& \\(n_{\\piter}\\) ");
fprintf("& \\(e_{\\piter}\\) ");
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

t_offline = result.t_construct_BF + result.t_construct_NS;
fprintf("\\(2^{%d}\\) ", log2(result.N));
fprintf("& %s ", form);
fprintf("& %.1e ", result.t_gmres);
fprintf("& %d ", result.num_iter_gmres);
fprintf("& %.1e ", result.rel_err_gmres);
fprintf("& %.1e ", t_offline);
fprintf("& %.1e ", result.t_pgmres);
fprintf("& %d ", result.num_iter_pgmres);
fprintf("& %.1e ", result.rel_err_pgmres);
fprintf("\\\\ \n");

end
