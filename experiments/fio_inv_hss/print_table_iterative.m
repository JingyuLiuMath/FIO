function print_table_iterative(...
    result_list, ...
    caption_name, label_name, ...
    tol_hss_display_list)

num_n = size(result_list, 1);
num_hss_tol = size(result_list, 2);

fprintf("\n");

fprintf("\\begin{table}[tbhp]\n");
fprintf("\\centering\n");
% fprintf("\\begin{tabular}{c | ccc | cccc}\n");
fprintf("\\begin{tabular}{c | ccc");
for it_hss_tol = 1 : num_hss_tol
    fprintf(" | cccc");
end
fprintf("}\n");
fprintf("\\toprule\n");

fprintf("\\multirow{2}{*}{\\(\\numtot\\)} ");
fprintf("& \\multirow{2}{*}{\\(t_{\\iter}\\)} ");
fprintf("& \\multirow{2}{*}{\\(n_{\\iter}\\)} ");
fprintf("& \\multirow{2}{*}{\\(e_{\\iter}\\)} ");

for it_hss_tol = 1 : num_hss_tol
    tol_hss_display = tol_hss_display_list(it_hss_tol);
    fprintf("& \\multicolumn{4}{c}{\\(\\varepsilon = %s\\)}", tol_hss_display);
    fprintf("\\\\ \n");
end

for it_hss_tol = 1 : num_hss_tol
    fprintf("& & & & ");
    fprintf("\\(t_{\\pre}\\) ");
    fprintf("& \\(t_{\\piter}\\) ");
    fprintf("& \\(n_{\\piter}\\) ");
    fprintf("& \\(e_{\\piter}\\) ");
    fprintf("\\\\ \n");
end

for it_n = 1 : num_n
    curr_result = result_list(it_n, :);
    print_result_iterative(curr_result, tol_hss_display_list);
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\end{table}\n");
fprintf("\n");


end

function print_result_iterative(result_list, tol_hss_display_list)

num_hss_tol = size(result_list, 2);
fprintf("\\midrule\n");
if isfield(result_list(1), "n")
    n = result_list(1).n;
    tmp_str = string(n) + "^{2}";
else
    N = result_list(1).N;
    p = log2(N);
    tmp_str = "2^{" + string(p) + "}";
end
fprintf("\\(%s\\) ", tmp_str);

fprintf("& %.1e ", result_list(1).t_cg);
fprintf("& %d ", result_list(1).iter_cg);
fprintf("& %.1e ", result_list(1).rel_err_cg);
for it_hss_tol = 1 : num_hss_tol
    t_pre = result_list(it_hss_tol).t_construct_BF ...
        + result_list(it_hss_tol).t_construct_HSS ...
        + result_list(it_hss_tol).t_factor_HSS;
    fprintf("& %.1e ", t_pre);
    fprintf("& %.1e ", result_list(it_hss_tol).t_pcg);
    fprintf("& %d ", result_list(it_hss_tol).iter_pcg);
    fprintf("& %.1e ", result_list(it_hss_tol).rel_err_pcg);
    fprintf("\\\\ \n");
end

end