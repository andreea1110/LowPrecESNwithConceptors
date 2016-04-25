function res = set_precision_sc(parameters)
%%% parameters.D = data matrix
%%% parameters.b = number of admitted values for the elements of D
%%% e.g: parameters.b = 3 => the elements in D can be -c, 0 or c, for a
%%% real constant c;
%%%% Return: D_lp = data matrix with precision log2(2*parameters.b + 1)

D = parameters.D;
noClusters = floor(parameters.b/2);
lowPrecD = zeros(size(D));
Dvec = abs(reshape(D', 1, numel(D)));
sigma = 1;
W = SimGraph_Full(Dvec, sigma);
[C, L, U] = SpectralClustering(W, noClusters, 1);

for i = 1:length(Dvec)
    A{C(i)} = Dvec(i);
end

centers = zeros(1, noClusters);

for i = 1:noClusters
    centers(i) = mean(A{i});
end

for i = 1:size(D, 1)
    for j = 1:size(D, 2)
        dif = zeros(1, noClusters);
        for k = 1:noClusters
            dif(k) = abs(abs(D(i, j)) - centers(k));
        end
        [valMin, idx2] = min(dif);
        lowPrecD(i, j) = centers(idx2);
    end
end

res.D_lp = sign(D).*lowPrecD;
end