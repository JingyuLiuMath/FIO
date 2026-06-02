function y = apply_mbf_batch(Factors, y)

N = size(y, 1);
M = 1;
for i = 1:size(Factors,1)-1
    M = max(M, size(Factors{i,1}.M, 1));
end

num_cols = size(y, 2);

memory_limit_gb = 128;
memory_limit_complex_double = memory_limit_gb * 1024^3 / 8 / 2;

num_cols_limit = floor(memory_limit_complex_double/ M);

batch_size = floor(num_cols / num_cols_limit);
res_size = num_cols - batch_size * num_cols_limit;

batch_size_list = num_cols_limit * ones(batch_size, 1);

if res_size ~= 0
    batch_size_list = [...
        batch_size_list;
        res_size];
end

num_batch = size(batch_size_list, 1);

if num_batch > 1
    fprintf("      num of cols: %d\n", num_cols);
    fprintf("      num_cols_limit: %d\n", num_cols_limit);
    fprintf("      batch_size: %d\n", batch_size);
end

offset = 0;
for it = 1 : num_batch
    curr_batch_size = batch_size_list(it);
    curr_ind = (offset + 1) : (offset + curr_batch_size);
    y(:, curr_ind) = apply_mbf(Factors, y(:, curr_ind));
    offset = offset + curr_batch_size;
end

end