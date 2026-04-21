function s = ShowSingularValueDecay(A, title_name)
% ShowSingularValueDecay

% Jingyu Liu, October 23, 2023.

arguments (Input)
    A (:, :) double;
    title_name string = "ShowSingularValueDecay";
end

arguments (Output)
    s (:, 1) double;
end

s = svd(A);
figure();
semilogy(1 : length(s), s / s(1), "LineWidth", 2);
xlabel("i");
ylabel("\sigma_i / \sigma_1");
ylim([1e-10, 1]);
title(title_name);

end