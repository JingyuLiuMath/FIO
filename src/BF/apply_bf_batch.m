function y = apply_bf_batch(Factor, y)

M = max(1, size(Factor.SigmaM, 1));
num_cols = size(y, 2);

memory_limit_gb = 128;
memory_limit_complex_double = ...
    memory_limit_gb * 1024^3 / 8 / 2;
num_cols_limit = max(...
    1, floor(memory_limit_complex_double / M));

num_full_batches = floor(num_cols / num_cols_limit);
num_remaining_cols = num_cols ...
    - num_full_batches * num_cols_limit;
batch_size_list = num_cols_limit ...
    * ones(num_full_batches, 1);
if num_remaining_cols ~= 0
    batch_size_list = [batch_size_list; num_remaining_cols];
end

offset = 0;
for it = 1 : length(batch_size_list)
    curr_batch_size = batch_size_list(it);
    curr_ind = offset + (1 : curr_batch_size);
    y(:, curr_ind) = ApplyFactor(Factor, y(:, curr_ind));
    offset = offset + curr_batch_size;
end

function y = ApplyFactor(Factor, y)

y = Factor.V' * y;
for i = length(Factor.BTol) : -1 : 1
    y = Factor.BTol{i}' * y;
end
y = Factor.SigmaM * y;
for i = 1 : length(Factor.ATol)
    y = Factor.ATol{i} * y;
end
y = Factor.U * y;

end

end
