function y = apply_fbf_adj(Factor, y)

y = Factor.U'*y;

for i=length(Factor.GTol):-1:1
    y = Factor.GTol{i}'*y;
end

y = Factor.M'*y;

for i=1:length(Factor.HTol)
    y = Factor.HTol{i}'*y;
end

y = Factor.V'*y;

end
