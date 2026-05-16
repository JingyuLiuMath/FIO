function [Omega, Y] = BBC_Indep_Y_Omega(A, s)

arguments (Input)
    A HSS_Herm;
    s (1, 1) double;
end

arguments (Output)
    Omega (:, :) double;
    Y (:, :) double;
end

if A.leaf_ == 1
    Omega = A.BBC_Omega_;
    Y = A.BBC_Y_;
    A.BBC_Omega_ = [];
    A.BBC_Y_ = [];
else
    Omega = zeros(A.level_size_, s);
    Y = zeros(A.level_size_, s);
    offset_Om = 0;
    offset_Y = 0;
    for i = 1 : A.num_children_
        curr_size_Om = size(A.children_{i}.BBC_Omega_, 1);
        Omega((offset_Om + 1) : (offset_Om + curr_size_Om), : )...
            = A.children_{i}.BBC_Omega_;
        A.children_{i}.BBC_Omega_ = [];
        offset_Om = offset_Om + curr_size_Om;

        curr_size_Y = size(A.children_{i}.BBC_Y_, 1);
        Y((offset_Y + 1) : (offset_Y + curr_size_Y), :) ...
            = A.children_{i}.BBC_Y_;
        A.children_{i}.BBC_Y_ = [];
        offset_Y = offset_Y + curr_size_Y;
    end
end

end