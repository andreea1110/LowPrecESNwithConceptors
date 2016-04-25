function res = set_precision_pca(parameters)
%%% parameters.D = data matrix 
%%% parameters.b = number of admitted values for the elements of D
%%% Return: D_lp = data matrix with precision log2(2*parameters.b + 1), 
%%%         maxd = the resulting discretization element
%%%         uv = number of unique values in D_lp

D = parameters.D;
nrows = size(D, 1);
ncols = size(D, 2);
D_save = D;

if nrows == 1 || ncols == 1
    if nrows == 1 
        D = reshape(D, sqrt(ncols), sqrt(ncols));
    else
        D = reshape(D, sqrt(nrows), sqrt(nrows));
    end
end
nrows = size(D, 1);
ncols = size(D, 2);

normD = (D - repmat(mean(D, 1), [nrows, 1]));
corrD = zeros(ncols, ncols);

for i = 1:ncols
    for j = i:ncols
        corrD(i, j) = sum(normD(:, i).*normD(:, j))/nrows;
        corrD(j, i) = sum(normD(:, i).*normD(:, j))/nrows;
    end
end
[eigvec, eigval] = eig(corrD);
[val, idx] = sort(diag(eigval), 'descend');
peigvec = eigvec(:, idx(1));
projD = peigvec'*normD;
ua_projD = unique(abs(projD));

b = parameters.b;
maxd = max(ua_projD);
D_lp = sign(D_save) .* ((round(abs(D_save) * b / maxd))*(maxd / b));
uv = length(unique(D_lp));

res.D_lp = D_lp;
res.maxd = maxd;
res.uv = uv;
end