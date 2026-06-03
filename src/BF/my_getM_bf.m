function M = my_getM_bf(BF, type_bf)
switch type_bf
    case "bf"
        M = size(BF.M, 1);
    case "mbf"
        M = 1;
        for i = 1:size(BF,1)-1
            M = max(M, size(BF{i,1}.M, 1));
        end
end

end