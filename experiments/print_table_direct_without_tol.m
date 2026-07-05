function print_table_direct_without_tol(result_list, ...
    caption_name, label_name)

num_n = size(result_list, 1);

fprintf("\n");

fprintf("\\begin{table}[tbhp]\n")
fprintf("\\centering\n")
fprintf("\\begin{tabular}{c | cc | cc | cc | c}\n")
fprintf("\\toprule\n")

fprintf("\\(\\numtot\\) ");

fprintf("& \\(t_{\\construct \\BF}\\) ");
fprintf("& \\(e_{\\construct \\BF}\\) ");

fprintf("& \\(t_{\\construct \\HSS}\\) ");
fprintf("& \\(e_{\\construct \\HSS}\\) ");

fprintf("& \\(t_{\\factor \\HSS}\\) ");
fprintf("& \\(t_{\\solve \\HSS}\\) ");
fprintf("& \\(e_{\\solve}\\) ");

fprintf("\\\\ \n");
for it_n = 1 : num_n
    curr_result = result_list(it_n, :);
    print_result_direct(curr_result);
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\end{table}\n");
fprintf("\n");


end

function print_result_direct(result_list)

num_hss_tol = 1;
fprintf("\\midrule\n");
if isfield(result_list, "n")
    n = result_list.n;
    tmp_str = string(n) + "^{2}";
else
    N = result_list.N;
    p = log2(N);
    tmp_str = "2^{" + string(p) + "}";
end
fprintf("\\(%s\\) ", tmp_str);


for it_hss_tol = 1 : num_hss_tol    
    fprintf("& %.1e ", result_list.t_construct_BF);
    fprintf("& %.1e ", result_list.rel_err_BF);
    
    fprintf("& %.1e ", result_list.t_construct_HSS);
    fprintf("& %.1e ", result_list.rel_err_HSS);
    
    fprintf("& %.1e ", result_list.t_factor_HSS);
    fprintf("& %.1e ", result_list.t_direct);
    fprintf("& %.1e ", result_list.rel_err_direct);
    
    fprintf("\\\\ \n");
end

end