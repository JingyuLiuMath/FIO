function U = Matrix_U(A)
% Matrix_U

arguments (Input)
    A HSS_Herm;
end

arguments (Output)
    U (:, :) double;
end

if A.leaf_ == 1
    U = A.Umat_;
else
    U = zeros(A.size_, A.rank_);
    offset = 0;
    for i = 1 : A.num_children_
        current_size = A.children_{i}.size_;
        U((offset + 1) : (offset + current_size), :) ...
            = A.children_{i}.Matrix_U() * A.Rmat_{i};
        offset = offset + current_size;
    end
end


end