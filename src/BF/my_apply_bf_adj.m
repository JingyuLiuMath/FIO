function y = my_apply_bf_adj(BF, y, type_bf)

switch type_bf
    case "bf"
        y = apply_fbf_adj_batch(BF, y);
    case "mbf"
        y = apply_mbf_adj_batch(BF, y);
end

end
