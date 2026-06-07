function y = my_apply_bf(BF, y, type_bf)

switch type_bf
    case "bf"
        y = apply_fbf_batch(BF, y);
    case "mbf"
        y = apply_mbf_batch(BF, y);
    case "pbf"
        y = apply_fbf_batch(BF, y);
    case "mybf"
        y = BF.Apply(y);
end

end
