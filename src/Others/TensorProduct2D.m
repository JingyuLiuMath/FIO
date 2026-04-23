function x_pos = TensorProduct2D(x1_range, x2_range)
% TensorProduct2D

% Jingyu Liu, January 27, 2025.

arguments (Input)
    x1_range (:, 1) double;
    x2_range (:, 1) double;
end

arguments (Output)
    x_pos (:, 2) double;
end

[x1_pos, x2_pos] = ndgrid(...
    x1_range, x2_range);
x1_pos = x1_pos(:);
x2_pos = x2_pos(:);
x_pos = [x1_pos, x2_pos];

end