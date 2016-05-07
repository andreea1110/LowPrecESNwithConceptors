function res = set_precision_ga(parameters)
%%% Reduce the precision of matrix D to 3 values -c, 0, c, where
%%% c is determined by using a genetic algorithm
%%% Return: D_lp = data matrix with precision log2(3)
D = parameters.D;
D_abs = abs(parameters.D);
nvars = 2;    % Number of variables
LB = [0 0];   % Lower bound
UB = [0.4 0.2];  % Upper bound -> UB(1) = max_psi; UB(2) = max_phi
ObjectiveFunction = @(x)fitness_function(x, parameters);
options = gaoptimset('Display', 'iter', 'PopulationSize', 20, 'Generations', 100, 'MigrationInterval', 20);
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[], options)

psi_opt = x(1);
phi_opt = x(2);
D_abs(D_abs >= phi_opt) = psi_opt;
D_abs(D_abs < phi_opt) = 0;
res.D_lp = sign(D).*D_abs;
end