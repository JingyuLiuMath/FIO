function y = apply_mbf_adj(Factors, x)

if( size(x,2)==0 )
    return;
end

y = zeros(size(x));

for i=1:size(Factors,1)-1
    yi = apply_fbf_adj(Factors{i,1},x);
    y(Factors{i,2}, :) = yi;
end
yi = Factors{end,1}' * x;
y(Factors{end,2},:) = yi;

end