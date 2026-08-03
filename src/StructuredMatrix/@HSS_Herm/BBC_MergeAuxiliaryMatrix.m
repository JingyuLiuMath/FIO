function BBC_MergeAuxiliaryMatrix(A)
% BBC_MergeAuxiliaryMatrix

arguments (Input)
    A HSS_Herm;
end

A.BBC_Omega_ = [];
A.BBC_Y_ = [];
for i = 1 : A.num_children_
    A.BBC_Omega_ = [A.BBC_Omega_;...
        A.children_{i}.Umat_' * A.children_{i}.BBC_Omega_];
    A.BBC_Y_ = [A.BBC_Y_; ...
        A.children_{i}.Umat_' ...
        * (A.children_{i}.BBC_Y_ ...
        - A.children_{i}.BBC_A_ * A.children_{i}.BBC_Omega_)];

    % Clear.
    A.children_{i}.BBC_Omega_ = [];
    A.children_{i}.BBC_Y_ = [];
    A.children_{i}.BBC_A_ = [];
    if A.children_{i}.leaf_ == 0
        A.children_{i}.Umat_ = [];
    end

end
