function M = my_getM_bf(BF, type_bf)
switch type_bf
    case "bf"
        M = size(BF.SigmaM, 1);
    case "fbf"
        M = size(BF.M, 1);
    case "mbf"
        M = 1;
        for i = 1 : (size(BF, 1) - 1)
            M = max(M, size(BF{i, 1}.SigmaM, 1));
        end
    case "pbf"
        M = size(BF.SigmaM, 1);
    case "fmbf"
        M = 1;
        for i = 1 : (size(BF, 1) - 1)
            M = max(M, size(BF{i, 1}.M, 1));
        end
    case "fpbf"
        M = size(BF.M, 1);
    case "mybf"
        M = BF.N_;
    otherwise
        error("FIO:BFGetM:UnsupportedType", ...
            "Unsupported BF type: %s.", type_bf);
end

end
