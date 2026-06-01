function print_table_direct(result_list, ...
    caption_name, label_name, ...
    tol_hss_display_list)

num_n = size(result_list, 1);

fprintf("\n");

fprintf("\\begin{table}[tbhp]\n")
fprintf("\\centering\n")
fprintf("\\begin{tabular}{c c | cc | cc | ccc}\n")
fprintf("\\toprule\n")

fprintf("\\(N\\) ");

fprintf("& \\(\\varepsilon\\) ");

fprintf("& \\(t_{\\mathrm{cBF}}\\) ");
fprintf("& \\(e_{\\mathrm{BF}}\\) ");

fprintf("& \\(t_{\\mathrm{cHSS}}\\) ");
fprintf("& \\(e_{\\mathrm{HSS}}\\) ");

fprintf("& \\(t_{\\mathrm{fHSS}}\\) ");
fprintf("& \\(t_{\\mathrm{sHSS}}\\) ");
fprintf("& \\(e_{\\mathrm{s}}\\) ");

fprintf("\\\\ \n");
for it_n = 1 : num_n
    curr_result = result_list(it_n, :);
    print_result_direct(curr_result, tol_hss_display_list);
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\end{table}\n");
fprintf("\n");


end

function print_result_direct(result_list, tol_hss_display_list)

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
fprintf("\\multirow{%d}{*}{\\(%s\\)} ", num_hss_tol, tmp_str);

for it_hss_tol = 1 : num_hss_tol
    tol_hss = result_list(it_hss_tol).tol_hss;
    tol_hss_display = tol_hss_display_list(it_hss_tol);
    if it_hss_tol == 1
        fprintf(" & \\(%s\\) ", tol_hss_display);
    else
        fprintf(" & & \\(%s\\) ", tol_hss_display);
    end
    
    fprintf("& %.1e ", result_list(it_hss_tol).t_construct_BF);
    fprintf("& %.1e ", result_list(it_hss_tol).rel_err_BF);
    
    fprintf("& %.1e ", result_list(it_hss_tol).t_construct_HSS);
    fprintf("& %.1e ", result_list(it_hss_tol).rel_err_HSS);
    
    fprintf("& %.1e ", result_list(it_hss_tol).t_factor_HSS);
    fprintf("& %.1e ", result_list(it_hss_tol).t_direct);
    fprintf("& %.1e ", result_list(it_hss_tol).rel_err_direct);
    
    fprintf("\\\\ \n");
end

end