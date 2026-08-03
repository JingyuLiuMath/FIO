function Construct(BF, a_func, phi_func, r, tol)
% Construct constructs and compresses an interpolation BHP.

% Jingyu Liu, August 3, 2026.

arguments (Input)
    BF Butterfly;
    a_func function_handle;
    phi_func function_handle;
    r (1, 1) double;
    tol (1, 1) double;
end

BF.constructed_ = false;
BF.Construct_Interpolation(a_func, phi_func, r);
BF.OutCompression(tol);
BF.InCompression(tol);
BF.constructed_ = true;

end
