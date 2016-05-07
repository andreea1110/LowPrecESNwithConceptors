function y = fitness_function(x, parameters)
%%% Function which we want to minimize by the genetic algorithm
%%% x(1) = psi; x(2) = phi
D = parameters.D;
D_abs = abs(D);
psi = x(1);
phi = x(2);
D_abs(D_abs >= phi) = psi;
D_abs(D_abs < phi) = 0;
D = sign(D).*D_abs;

if strcmp(parameters.mat, 'W')
    parameters.ESN.W = D;
elseif strcmp(parameters.mat, 'Wout')
    parameters.ESN.Wout = D;
elseif strcmp(parameters.mat, 'C1')
    parameters.C.mat{1} = D;
elseif strcmp(parameters.mat, 'C2')
    parameters.C.mat{2} = D;
elseif strcmp(parameters.mat, 'C3')
    parameters.C.mat{3} = D;
elseif strcmp(parameters.mat, 'C4')
    parameters.C.mat{4} = D;
else
    
end

netOutput = testESNC(parameters.ESN, parameters.C, parameters.patterns);
[netOutputInt, patternsInt] = interpolate(netOutput, parameters.train_pPL, parameters.intRate, parameters.ESN);
[netOutputIntAl, patternsIntResized] = align_signals(netOutputInt, patternsInt, parameters.intRate, parameters.ESN);
y = compute_error(netOutputIntAl, patternsIntResized, size(patternsIntResized, 2));
end