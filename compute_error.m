function meanNRMSE = compute_error(netOutputInt, patternsInt, len)
% Computes and prints the normalized root mean squared error (NRMSE)
% between the output signal and the target signal, for each
% of the signals; also computes the mean error for all the signals
% (patterns)
    k = size(patternsInt, 1);
    NRMSE_vec = zeros(1, k);
    %disp('**********************');
    for i = 1:k
        NRMSE_i = nrmse(netOutputInt(i, 1:len), patternsInt(i, 1:len));
        NRMSE_vec(i) = NRMSE_i;    
        %fprintf('NRMSE pattern %d: %g\n', i, NRMSE_i);
    end
    %fprintf('mean NRMSE: %g\n', mean(NRMSE_vec));
    meanNRMSE = mean(NRMSE_vec);
end