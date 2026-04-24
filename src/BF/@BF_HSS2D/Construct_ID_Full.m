function Construct_ID_Full(G, exp_phi_func, x, xi, rank_or_tol)
% Construct

arguments (Input)
    G BF_HSS2D;
    exp_phi_func function_handle;
    x (:, 2) double;
    xi (:, 2) double;
    rank_or_tol (1, 1) double;
end

inact = [];
for level = G.max_level_ : -1 : 1
    G.ConstructGenerators_ID_Full(...
        exp_phi_func, ...
        x, xi, ...
        inact, ...
        level, rank_or_tol);
    inact = G.Deactivate(level, inact);
end

G.ConstructRootGenerators2(exp_phi_func, x, xi);

end