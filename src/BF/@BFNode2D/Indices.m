function ind = Indices(B)

arguments (Input)
    B BFNode2D;
end

arguments (Output)
    ind (:, 1) double;
end

ind1 = (B.node1_.ind_start_ : B.node1_.ind_end_)' + 1;
ind2 = (B.node2_.ind_start_ : B.node2_.ind_end_)' + 1;
ind = TensorProduct2D(ind1, ind2);
ind = sub2ind([B.n_, B.n_], ind(:, 1), ind(:, 2));

end