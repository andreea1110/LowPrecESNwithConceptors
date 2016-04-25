function netOutput = testESNC(ESN, C, patterns)
 % Computes the ouput patterns, using the trained ESN with conceptors
    k = length(patterns);
    netOutput = zeros(k, ESN.testLength);   
    for i = 1:k
        Ci = C.mat{i}; % store conceptor matrix for pattern i
        x = 0.5 * randn(ESN.DR.size,1); % initialize the activations of the units in the DR
        for j = 1:ESN.testLength
            x = Ci * tanh(ESN.W *  x + ESN.Wbias)+ ESN.stateNL * (rand(ESN.DR.size,1)-0.5);
            netOutput(i, j) = ESN.Wout * x;
        end     
    end
end