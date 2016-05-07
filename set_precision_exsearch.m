function res = set_precision_exsearch(parameters)
%%% Reduce the precision of matrix D to 3 values -c, 0, c, where
%%% c is determined by exhaustive search
%%% Return: D_lp = data matrix with precision log2(3)
D = parameters.D;
mat = parameters.mat;

%phi_range = unique(abs(D));
%psi_range = unique(abs(D));
nvals = 45; % number of entries from original D to do exhaustive search on

phi_range = linspace(0, 0.2, nvals);
psi_range = linspace(0, 0.4, nvals);
% idx = [];
% while length(idx) < nvals
%     idx = [idx, randint(1, 1, [1, length(unique(abs(D)))])];
%     idx = unique(idx);
% end
%phi_range = sort(phi_range(idx));
%psi_range = sort(psi_range(idx));

error = zeros(nvals, nvals);

min_err = Inf;
D_save = D;
D_abs = abs(D);
phi_opt = 0;
psi_opt = 0;
est_run_time = 0;
for i = 1:length(phi_range)
    if parameters.verbose
        tic;
    end
    phi = phi_range(i);
    if parameters.verbose
        fprintf('%d) ', i)
        fprintf('phi = %0.3g\n', phi)
        fprintf('min_err = %0.3g\n', min_err)
    end
    for j = 1:length(psi_range)
        psi = psi_range(j);
        
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
        error(i, j) = compute_error(netOutputIntAl, patternsIntResized, size(patternsIntResized, 2));
        if error(i, j) < min_err
            min_err = error(i, j);
            phi_opt = phi;
            psi_opt = psi;
        end
        D = D_save;
        D_abs = abs(D);
    end
    
    if parameters.verbose
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
end

save(strcat('variables/params', mat, '.mat'), 'psi_range', 'phi_range', 'error');

if parameters.verbose
    fprintf('phi_opt (%s): %d\n', mat, phi_opt);
    fprintf('psi_opt (%s): %d\n', mat, psi_opt);
    fprintf('mean of vals in original %s: %d\n', mat, mean(D(:)));
    fprintf('std of vals in original %s: %d\n', mat, std(D(:)));
    
    % plot distributions of original matrices, together with the optimal
    % psi and phi
    f = figure(); clf;
    set(gca,'fontsize',16);
    hold on;
    hist(D(:), 40)
    ylim = get(gca,'ylim');
    h1 = line([phi_opt phi_opt], ylim, 'Color','g', 'Linewidth', 2);
    line([-phi_opt -phi_opt], ylim, 'Color','g', 'Linewidth', 2);
    h2 = line([psi_opt psi_opt], ylim, 'Color','r', 'Linewidth', 2);
    line([-psi_opt -psi_opt], ylim, 'Color','r', 'Linewidth', 2);
    legend([h1, h2], {strcat('|\phi| = ', sprintf('%g', phi_opt)), strcat('|\psi| = ', sprintf('%g', psi_opt))});
    title(strcat(mat,sprintf('\n mean = %g, var = %g, std = %g', mean(D(:)), var(D(:)), std(D(:)))))
    xlabel('values');
    ylabel('counts');
    saveas(f, strcat('images/', mat, '.png'), 'png');
    saveas(f, strcat('images/', mat, '.fig'), 'fig');
    
    % plot error vs. psi and phi
    f = figure();clf; 
    surf(psi_range, phi_range, error);
    colormap(jet);
    xlabel('\psi', 'fontsize',16);
    ylabel('\phi', 'fontsize',16);
    zlabel('NRMSE', 'fontsize',16);
    title(['Error vs. optimization parameters for ', mat], 'fontsize',16);
    hold on;
    set(gca,'fontsize',16);
    grid on;
    saveas(f, strcat('images/error_landscape_', mat, '.png'), 'png');
    saveas(f, strcat('images/error_landscape_', mat, '.fig'), 'fig');
    
    % plot error vs. psi and phi - contour plot
    f = figure();clf; 
    contourf(psi_range, phi_range, error);
    colormap(jet);
    %clabel(c,h);  
    xlabel('\psi', 'fontsize',16);
    ylabel('\phi', 'fontsize',16);
    %zlabel('NRMSE', 'fontsize',16);
    title(['Error contour plot for ', mat], 'fontsize',16);
    hold on;
    set(gca,'fontsize',16);
    grid on;
    saveas(f, strcat('images/error_contour', mat, '.png'), 'png');
    saveas(f, strcat('images/error_contour', mat, '.fig'), 'fig');
    
end

D_abs(D_abs >= phi_opt) = psi_opt;
D_abs(D_abs < phi_opt) = 0;
res.D_lp = sign(D).*D_abs;

if parameters.verbose
    fprintf('no unique vals in %s: %d\n', mat, length(unique(res.D_lp)));
    fprintf('unique vals in %s: %s\n', mat, sprintf('%d ', unique(res.D_lp)));
end
end