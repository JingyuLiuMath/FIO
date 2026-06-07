function res = my_nnz_bf(BF, type_bf, verbose)

arguments (Input)
    BF;
    type_bf string;
    verbose (1, 1) double = 0;
end

switch type_bf
    case "bf"
        res = my_nnz_bf_local(BF);
        if verbose == 1
            print_nnz_local(res);
        end
    case "mbf"
        res = cell(size(BF,1), 1);
        for i = 1 : (size(BF, 1) - 1)
            res{i} = my_nnz_bf_local(BF{i, 1});
        end
        res{size(BF, 1)} = nnz(BF{size(BF, 1), 1});
        if verbose == 1
            print_nnz_mbf(res);
        end
end

end

function res = my_nnz_bf_local(BF)

res = struct();

res.nnz_U = nnz(BF.U);

res.nnz_GTol = zeros(size(BF.GTol));
for i = 1 : size(BF.GTol)
    res.nnz_GTol(i) = nnz(BF.GTol{i});
end

res.nnz_M = nnz(BF.M);

res.nnz_HTol = zeros(size(BF.HTol));
for i = 1 : size(BF.HTol)
    res.nnz_HTol(i) = nnz(BF.HTol{i});
end

res.nnz_V = nnz(BF.V);

end

function print_nnz_local(nnz_BF)

fprintf("  nnz_U: %d\n", nnz_BF.nnz_U);
for i = 1 : size(nnz_BF.nnz_GTol)
    fprintf("  nnz_GTol(%d): %d\n", i, nnz_BF.nnz_GTol(i));
end
fprintf("  nnz_M: %d\n", nnz_BF.nnz_M);
for i = 1 : size(nnz_BF.nnz_HTol)
    fprintf("  nnz_HTol(%d): %d\n", i, nnz_BF.nnz_HTol(i));
end
fprintf("  nnz_V: %d\n", nnz_BF.nnz_V);

end

function print_nnz_mbf(nnz_mbf)

m = size(nnz_mbf, 1);
for i = 1 : (m - 1)
    fprintf("the %d-th row in mbf\n", i);
    print_nnz_local(nnz_mbf{i});
end
i = m;
fprintf("the %d-th row in mbf\n", i);
fprintf("  nnz: %d\n", nnz_mbf{i});

end