function res = set_precision_distr1(parameters)
%%% Reduce the precision of the matrix/vector D = (d_ij).
%%% parameters.D = data matrix
%%% parameters.b = number of admitted values for the elements of D
%%% e.g: b = 3 => the elements in D can be -c, 0 or c, for a
%%% real constant c;
%%% Return: D_lp = data matrix with precision log2(2*parameters.b + 1)

b = parameters.b;
D = parameters.D;
lowPrecD = zeros(size(D)); % initializing the low precision matrix
[counts, centers] = hist(abs(D(:)), b); % histogram with b bins; counts = number of elements per bin; centers = centers of the bins

% iterate through the D matrix
for i = 1:size(D, 1)
    for j = 1:size(D, 2)
        dif = zeros(1, b); % vector of the differences between entry d_ij and each of the centers
        for k = 1:b
            dif(k) = abs(abs(D(i, j)) - centers(k));
        end
        [valMin, idx2] = min(dif); % setting d_ij to the center of the bin it fell into
        lowPrecD(i, j) = centers(idx2);
    end
end
res.D_lp = sign(D).*lowPrecD; % reassigning signs
end