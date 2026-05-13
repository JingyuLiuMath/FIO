function S = sparse_blkdiag(varargin)

n_blocks = nargin;
rows = cell(n_blocks, 1);
cols = cell(n_blocks, 1);
vals = cell(n_blocks, 1);

row_offset = 0;
col_offset = 0;

for k = 1:n_blocks
    B = varargin{k};
    [m, n] = size(B);

    [i, j, v] = find(B);

    rows{k} = i + row_offset;
    cols{k} = j + col_offset;
    vals{k} = v;

    row_offset = row_offset + m;
    col_offset = col_offset + n;
end

S = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), ...
    row_offset, col_offset);
end