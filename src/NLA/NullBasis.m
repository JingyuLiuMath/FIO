function Z = NullBasis(B, k)
% NullBasis

% Jingyu Liu, December 4, 2024.

arguments (Input)
    B (:, :) double;
    k (1, 1) double
end

arguments (Output)
    Z (:, :) double;
end

[Z, ~, ~] = qr(B', "vector");
Z = Z(:, (end - k + 1) : 1 : end);

end