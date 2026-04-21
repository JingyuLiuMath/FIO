function BBC_FillAuxiliaryMatrix(A, Omega, Y)
% BBC_FillAuxiliaryMatrix

arguments (Input)
    A HSS_Herm;
    Omega (:, :) double;
    Y (:, :) double;
end

if A.leaf_ == 1
    A.BBC_Omega_ = Omega;
    A.BBC_Y_ = Y;
else
    % Recursion.
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.size_;
        A.children_{i}.BBC_FillAuxiliaryMatrix(...
            Omega((offset + 1) : (offset + current_size), :), ...
            Y((offset + 1) : (offset + current_size), :));
        offset = offset + current_size;
    end
end

end