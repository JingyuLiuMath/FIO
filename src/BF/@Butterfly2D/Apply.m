function f = Apply(BF, f)

arguments (Input)
    BF Butterfly2D;
    f (: ,:) double;
end

arguments (Output)
    f (:, :) double;
end

f = f(BF.xi_perm_, :);
f = bmatrix_mult_vec(BF.V_, f);

for ind_H = length(BF.H_) : -1 : 1
    f = reshape_H(BF.H_{ind_H}, f);
    f = bmatrix_mult_bvec(BF.H_{ind_H}, f);
end

f = reshape_H(BF.M_, f);
f = bmatrix_mult_bvec(BF.M_, f);
f = reshape_G(BF.M_, f);

for ind_G = length(BF.G_) : -1 : 1
    f = reshape_G(BF.G_{ind_G}, f);
    f = bmatrix_mult_bvec(BF.G_{ind_G}, f);
end

f = bmatrix_mult_bvec(BF.U_, f);
f = bvec2vec(f);
f = f(BF.x_perm_inv_, :);

end

function g = reshape_H(A, f)

arguments (Input)
    A (:, :) cell;
    f (:, :) cell;
end

arguments (Output)
    g (:, :) cell;
end

[mH, nH] = size(A);
[m, n] = size(f);

if mH == m && nH == n
    g = f;
    return;
end

g = cell(mH, nH);
for i_par = 1 : m
    for j = 1 : nH
        i_offset = 0;
        for i = [4 * i_par - 3, 4 * i_par - 2, 4 * i_par - 1, 4 * i_par]
            i_size = size(A{i, j}, 2);
            i_ind = (i_offset + 1) : (i_offset + i_size);
            y = 0;
            for j_ch = [4 * j - 3, 4 * j - 2, 4 * j - 1, 4 * j]
                y = y + f{i_par, j_ch}(i_ind, :);
            end
            i_offset = i_offset + i_size;
            g{i, j} = y;
        end
    end
end

end

function g = reshape_G(A, f)

arguments (Input)
    A (:, :) cell;
    f (:, :) cell;
end

arguments (Output)
    g (:, :) cell;
end

[mH, nH] = size(A);
[m, n] = size(f);

if mH == m && nH == n
    g = f;
    return;
end

g = cell(mH, nH);
for i_par = 1 : m
    for j = 1 : nH
        for i = [4 * i_par - 3, 4 * i_par - 2, 4 * i_par - 1, 4 * i_par]
            y = [];
            for j_ch = [4 * j - 3, 4 * j - 2, 4 * j - 1, 4 * j]
                y = [y; f{i_par, j_ch}];
            end
            g{i, j} = y;
        end
    end
end

end