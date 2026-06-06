function SetMidSize(BF)

arguments (Input)
    BF Butterfly2D;
end

mid_row_size = 0;
for i = 1 : size(BF.M_, 1)
    mid_row_size = mid_row_size + size(BF.M_{i, 1}, 1);
end

mid_col_size = 0;
for j = 1 : size(BF.M_, 2)
    mid_col_size = mid_col_size + size(BF.M_{1, j}, 2);
end

BF.mid_row_size_ = mid_row_size;
BF.mid_col_size_ = mid_col_size;

end