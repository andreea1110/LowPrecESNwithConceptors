function res = set_precision_rand(parameters)
%%% parameters.D = data matrix
%%% parameters.b = number of admitted values for the elements of D
%%% bound numerical precision of data D = (d_i) = parameters.D (number, vec, mat) to parameters.b linearly
%%% spaced values by:
%%% 1. choosing maxd randomly from the distribution obtained by the pca
%%% algorithm 'set_precision_pca'
%%% 2. splitting [0, maxd] into b equal intervals;
%%% 3. rounding all abs(d_i) toward the interval boundaries;
%%% 4. re-assigning signs
%%% The rounding is done on a linear data scale.
%%%% Return: D_lp = data matrix with precision log2(2*parameters.b + 1)
%%%          maxd = the resulting discretization element
%%%          uv = number of unique values in D_lp
b = parameters.b;
D = parameters.D;
[counts, centers] = hist(parameters.maxdvec, 30);
pmf = counts/sum(counts);
maxd = centers(discretesample(pmf, 1));
D_lp = sign(D) .* ((round(abs(D) * b / maxd))*(maxd / b));
uv = length(unique(D_lp));

res.D_lp = D_lp;
res.maxd = maxd;
res.uv = uv;
end
