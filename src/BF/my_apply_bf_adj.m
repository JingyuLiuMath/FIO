function y = my_apply_bf_adj(BF, y, type_bf)

switch type_bf
    case "bf"
        y = apply_bf_adj_batch(BF, y);
    case "fbf"
        y = apply_fbf_adj_batch(BF, y);
    case "mbf"
        y = apply_mbf_adj_batch(BF, y);
    case "pbf"
        y = apply_bf_adj_batch(BF, y);
    case "fmbf"
        y = apply_fmbf_adj_batch(BF, y);
    case "fpbf"
        y = apply_fbf_adj_batch(BF, y);
    case "mybf"
        y = BF.ApplyAdjoint(y);
    case "adjoint_view"
        y = conj(BF.alpha) * my_apply_bf(...
            BF.factor, y, BF.type_bf);
    otherwise
        error("FIO:BFApplyAdjoint:UnsupportedType", ...
            "Unsupported BF type: %s.", type_bf);
end

end
