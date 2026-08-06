function print_table_iterative(result_list, caption_name, label_name)

fprintf("\n");
fprintf("\\begin{table}[tbhp]\n");
fprintf("\\centering\n");
fprintf("\\begin{tabular}{c | ccc | cccc}\n");
fprintf("\\toprule\n");
fprintf("\\(\\numtot\\) ");
fprintf("& \\(t_{\\iter}\\) ");
fprintf("& \\(n_{\\iter}\\) ");
fprintf("& \\(e_{\\iter}\\) ");
fprintf("& \\(t_{\\offline}\\) ");
fprintf("& \\(t_{\\piter}\\) ");
fprintf("& \\(n_{\\piter}\\) ");
fprintf("& \\(e_{\\piter}\\) ");
fprintf("\\\\ \n");

for it_n = 1 : length(result_list)
    result = result_list{it_n};
    t_offline = result.t_construct_BF + result.t_construct_NS;
    fprintf("\\midrule\n");
    if isfield(result, "n")
        fprintf("\\(%d^{2}\\) ", result.n);
    else
        fprintf("\\(2^{%d}\\) ", log2(result.N));
    end
    fprintf("& %.1e ", result.t_gmres);
    fprintf("& %d ", result.num_iter_gmres);
    fprintf("& %.1e ", result.rel_err_gmres);
    fprintf("& %.1e ", t_offline);
    fprintf("& %.1e ", result.t_pgmres);
    fprintf("& %d ", result.num_iter_pgmres);
    fprintf("& %.1e ", result.rel_err_pgmres);
    fprintf("\\\\ \n");
end

fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\end{table}\n");
fprintf("\n");

end
