function err = HermError(A)

err = norm(A - A', "fro") / norm(A, "fro");

end