function RandInit(A, hss_rank_rule)

arguments (Input)
    A HSS_Herm;
    hss_rank_rule string;
end

for level = A.max_level_ : -1 : 1
    A.RandInit_Generators(level, hss_rank_rule);
end

A.RandInit_RootGenerators();

end