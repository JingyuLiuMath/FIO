function cnt = Nnz(BF)

arguments (Input)
    BF Butterfly;
end

arguments (Output)
    cnt (1, 1) double;
end

cnt = 0;
for i = 1 : size(BF.U_, 1)
    for j = 1 : size(BF.U_, 2)
        cnt = cnt + numel(BF.U_{i, j});
    end
end

for k = 1 : length(BF.G_)
    for i = 1 : size(BF.G_{k}, 1)
        for j = 1 : size(BF.G_{k}, 2)
            cnt = cnt + numel(BF.G_{k}{i, j});
        end
    end
end

for i = 1 : size(BF.M_, 1)
    for j = 1 : size(BF.M_, 2)
        cnt = cnt + numel(BF.M_{i, j});
    end
end

for k = 1 : length(BF.H_)
    for i = 1 : size(BF.H_{k}, 1)
        for j = 1 : size(BF.H_{k}, 2)
            cnt = cnt + numel(BF.H_{k}{i, j});
        end
    end
end

for i = 1 : size(BF.V_, 1)
    for j = 1 : size(BF.V_, 2)
        cnt = cnt + numel(BF.V_{i, j});
    end
end


end