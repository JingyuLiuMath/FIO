function mem = Storage(BF)

arguments (Input)
    BF Butterfly;
end

arguments (Output)
    mem (1, 1) double;
end

mem = 0;
for i = 1 : size(BF.U_, 1)
    for j = 1 : size(BF.U_, 2)
        mem = mem + byte_size(BF.U_{i, j});
    end
end

for k = 1 : length(BF.G_)
    for i = 1 : size(BF.G_{k}, 1)
        for j = 1 : size(BF.G_{k}, 2)
            mem = mem + byte_size(BF.G_{k}{i, j});
        end
    end
end

for i = 1 : size(BF.M_, 1)
    for j = 1 : size(BF.M_, 2)
        mem = mem + byte_size(BF.M_{i, j});
    end
end

for k = 1 : length(BF.H_)
    for i = 1 : size(BF.H_{k}, 1)
        for j = 1 : size(BF.H_{k}, 2)
            mem = mem + byte_size(BF.H_{k}{i, j});
        end
    end
end

for i = 1 : size(BF.V_, 1)
    for j = 1 : size(BF.V_, 2)
        mem = mem + byte_size(BF.V_{i, j});
    end
end


end