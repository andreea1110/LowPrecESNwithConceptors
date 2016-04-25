function [netOutputIntAl, patternsIntResized] = align_signals(netOutputInt, patternsInt, intRate, ESN)
%%% trying to allign the signals in phase as good as possible

    k = size(patternsInt, 1);
    L = size(netOutputInt, 2); % number of time steps in the output signal
    M = size(patternsInt, 2); % number of time steps in the input signal
    phasematches = zeros(k, L - M); 
    
    patternsIntResized = zeros(k, length(1:ESN.plotLength));
    netOutputIntAl = zeros(k, length(1:ESN.plotLength));
    
    for i = 1:k
        for phaseshift = 1:(L - M)
            phasematches(i, phaseshift) = norm(patternsInt(i, :) - netOutputInt(i, phaseshift:phaseshift+M-1));
        end
        [minVal minInd] = min(phasematches(i, :));
        netOutputIntAl(i, :) = netOutputInt(i, minInd:intRate:(minInd+intRate*ESN.plotLength-1));
        patternsIntResized(i, :) = patternsInt(i, 1:intRate:intRate*ESN.plotLength-1);   
    end
end