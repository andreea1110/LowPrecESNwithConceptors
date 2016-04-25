function [netOutputInt, patternsInt] = interpolate(netOutput, train_pPL, intRate, ESN)
    
    k = size(train_pPL, 2);
    patternsInt = zeros(k, length(1:(1/intRate):ESN.plotLength));
    netOutputInt = zeros(k, length(1:(1/intRate):ESN.testLength));
    
    for i = 1:k
        patternsInt(i, :) = interp1((1:ESN.plotLength)',train_pPL{1, i}', (1:(1/intRate):ESN.plotLength)', 'spline')';
        netOutputInt(i, :) = interp1((1:ESN.testLength)', netOutput(i, :)', (1:(1/intRate):ESN.testLength)', 'spline')';
    end

end
