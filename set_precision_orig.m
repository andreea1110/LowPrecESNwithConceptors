function res = set_precision_orig(parameters)
%%% bound numerical precision of data D = (d_i) = parameters.D (number, vec, mat) to parameters.b linearly
%%% spaced values by:
%%% 1. splitting [0, max(abs(D))] into b equal intervals;
%%% 2. rounding all abs(d_i) toward the interval boundaries;
%%% 3. re-assigning signs
%%% The rounding is done on a linear data scale.
%%% Return: D_lp = data matrix with precision log2(2*parameters.b + 1)

D = parameters.D;
b = parameters.b;
if b ~= Inf && max(max(abs(D))) > 0
    maxd = max(max(abs(D)));
    res.D_lp = sign(D) .* ((round(abs(D) * b / maxd))*(maxd / b));
end