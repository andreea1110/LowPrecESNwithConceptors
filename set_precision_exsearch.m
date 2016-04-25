function res = set_precision_exsearch(parameters)
%%% Reduce the precision of matrix D to 3 values -c, 0, c, where
%%% c is determined by exhaustive search
%%% Return: D_lp = data matrix with precision log2(3)

D = parameters.D;
mat = parameters.mat;
uniqD = unique(D);
fprintf('no unique el in W: %d\n', length(uniqD))
phi_range = unique(abs(D));
psi_range = unique(abs(D));
nvals = 100; % number of entries from original D to do exhaustive search on
idx = randint(1, nvals, [1, length(unique(abs(D)))]);
phi_range = phi_range(idx);
psi_range = psi_range(idx);

fprintf('length phi_range %d\n', length(phi_range));

err_vec = zeros(length(phi_range), length(psi_range));
min_err = Inf;
D_save = D;
D_abs = abs(D);
phi_opt = 0;
psi_opt = 0;
est_run_time = 0;
for i = 1:length(phi_range)
    tic;
    phi = phi_range(i);
    fprintf('%d) ', i)
    fprintf('phi = %0.3g\n', phi)
    fprintf('min_err = %0.3g\n', min_err)
    for j = 1:length(psi_range)
        psi = psi_range(j);
        
        D_abs(D_abs >= phi) = psi;
        D_abs(D_abs < phi) = 0;
        D = sign(D).*D_abs;
        %fprintf('no unique vals in W %d\n', length(unique(ESN.W)));
        
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
        err_vec(i, j) = compute_error(netOutputIntAl, patternsIntResized, size(patternsIntResized, 2));
        if err_vec(i, j) < min_err
            min_err = err_vec(i, j);
            phi_opt = phi;
            psi_opt = psi;
        end
        D = D_save;
        D_abs = abs(D);
    end
    
    % calculate and display run-time 
    elapsed_time = toc;
    if i == 1
        est_run_time = elapsed_time*length(phi_range);
        if est_run_time < 60
            fprintf('Estimated total run-time: %g seconds \n', est_run_time);
        else
            fprintf('Estimated total run-time: %g minutes \n', est_run_time/60);
        end
    end
    
    est_run_time = est_run_time - elapsed_time;
    if est_run_time < 60
        fprintf('Estimated remaining run-time: %g seconds \n', est_run_time);
    else
        fprintf('Estimated remaining run-time: %g minutes \n', est_run_time/60);
    end
end

fprintf('phi_opt (%s): %d\n', mat, phi_opt);
fprintf('psi_opt (%s): %d\n', mat, psi_opt);
fprintf('mean of vals in original %s: %d\n', mat, mean(D(:)));
fprintf('std of vals in original %s: %d\n', mat, std(D(:)));
f = figure();
hist(D(:), 40)
title(mat)
xlabel('values')
ylabel('counts')
saveas(f, strcat('images/', mat, '.png'), 'png')
saveas(f, strcat('images/', mat, '.fig'), 'fig')

D_abs(D_abs >= phi_opt) = psi_opt;
D_abs(D_abs < phi_opt) = 0;
res.D_lp = sign(D).*D_abs;
fprintf('no unique vals in %s: %d\n', mat, length(unique(lowPrecD)));
fprintf('unique vals in %s: %s\n', mat, sprintf('%d ', unique(lowPrecD)));

end