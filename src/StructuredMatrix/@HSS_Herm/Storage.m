function mem = Storage(A)

arguments (Input)
    A HSS_Herm;
end

arguments (Output)
    mem (1, 1) double;
end

mem = 0;

mem = mem + byte_size(A.Amat_);
mem = mem + byte_size(A.Umat_);

for i = 1:size(A.Rmat_, 2)
    mem = mem + byte_size(A.Rmat_{i});
end

for i = 1:size(A.Bmat_, 1)
    for j = 1:size(A.Bmat_, 2)
        mem = mem + byte_size(A.Bmat_{i, j});
    end
end

if A.leaf_ == 0
    for i = 1:A.num_children_
        mem = mem + A.children_{i}.Storage();
    end
end

end