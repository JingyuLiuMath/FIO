function r = Rank(A)
% Storage

arguments (Input)
    A HSS_Herm;
end

arguments (Output)
    r (1, 1) double;
end

r = A.rank_;

if A.leaf_ == 0
    for i = 1 : A.num_children_
        r = max(r, A.children_{i}.Rank());
    end
end

end