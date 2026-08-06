function num_nnz = my_nnz_bf(BF, type_bf)
% my_nnz_bf returns the total stored entries of a BF representation.

arguments (Input)
    BF;
    type_bf (1, 1) string;
end

arguments (Output)
    num_nnz (1, 1) double;
end

switch type_bf
    case {"bf", "pbf"}
        num_nnz = SmpNnz(BF);
    case {"fbf", "fpbf"}
        num_nnz = FbfNnz(BF);
    case "mbf"
        num_nnz = MultiscaleNnz(BF, @SmpNnz);
    case "fmbf"
        num_nnz = MultiscaleNnz(BF, @FbfNnz);
    case "mybf"
        num_nnz = BF.Nnz();
    case "adjoint_view"
        num_nnz = my_nnz_bf(BF.factor, BF.type_bf);
    otherwise
        error("FIO:BFNnz:UnsupportedType", ...
            "Unsupported BF type: %s.", type_bf);
end

end

function num_nnz = SmpNnz(Factor)

num_nnz = nnz(Factor.U) + nnz(Factor.SigmaM) ...
    + nnz(Factor.V);
for i = 1 : length(Factor.ATol)
    num_nnz = num_nnz + nnz(Factor.ATol{i});
end
for i = 1 : length(Factor.BTol)
    num_nnz = num_nnz + nnz(Factor.BTol{i});
end

end

function num_nnz = FbfNnz(Factor)

num_nnz = nnz(Factor.U) + nnz(Factor.M) + nnz(Factor.V);
for i = 1 : length(Factor.GTol)
    num_nnz = num_nnz + nnz(Factor.GTol{i});
end
for i = 1 : length(Factor.HTol)
    num_nnz = num_nnz + nnz(Factor.HTol{i});
end

end

function num_nnz = MultiscaleNnz(Factors, factor_nnz)

num_nnz = nnz(Factors{end, 1});
for i = 1 : (size(Factors, 1) - 1)
    num_nnz = num_nnz + factor_nnz(Factors{i, 1});
end

end
