function res = set_precision_distr2(parameters)
%%% Reduces the precision of matrix parameters.D to only
%%% 3 possible values: -c, 0, c
%%% Return: D_lp = data matrix with precision log2(3)

D = parameters.D;
stdevD = std(D(:));
meanD = mean(D(:));
D_abs = abs(D);

if strcmp(parameters.mat, 'W')
    phi = 0.8*stdevD;
    psi = 2.8*stdevD;
elseif strcmp(parameters.mat, 'Wout')
    phi = meanD;
    psi = 1.1*stdevD;
else
    phi = stdevD;
    %fprintf('max(D_abs(:))/stdevD = %g\n',max(D_abs(:))/stdevD)
    if max(D_abs(:))/stdevD > 5 
        psi = max(D_abs(:));
    else
        psi = 3*stdevD;
    end    
end

D_abs(D_abs >= phi) = psi;
D_abs(D_abs < phi) = 0;
res.D_lp = sign(D).*D_abs;

end