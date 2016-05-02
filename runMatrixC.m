function runMatrixC()
%%% Main runnable file.
%%% Set the desired method to reduce the bit precision of the parameters of
%%% the ESN with conceptors below, set the number of runs and adjust the
%%% plotting options and...press run. :-)
%%% Note: running the network at full precision is also an option.
%%% Possible ways to reduce the bit precision are:
%%% set_precision_orig  = 
%%% set_precision_distr1 = 
%%% set_precision_distr2 = 
%%% set_precision_pca = 
%%% set_precision_randunif =
%%% set_precision_exsearch = 
clear all;
clc;

%% Set options
settings.red_prec = true; % true = reduce the precision of the weight and conceptor matrices; false = run the algorithm at full precision
settings.red_prec_alg = 'set_precision_pca'; % choose the desired algorithm for reducing the precision
if settings.red_prec
    settings.b = 4; % the max number of unique values D_lp has; the actual max nr of values D_lp will have is 2*b + 1
else
    settings.b = Inf; % i.e. full precison 
end
settings.svd = false; % True = reduce precision of the SVD decomposition of all targeted matrices
                      % Note: The algorithms that work for the svd decomposition of the matrices
                      % are 'set_precision_orig' and 'set_precision_distr1'.
settings.plot_init_distr = false; % True = plot the initial distribution of all the weight and the conceptor matrices
settings.plot_signals = false; % True = plot the output of the neural network on top of the training signals;
settings.plot_full = false; % True = plot signals & 2 neurons & the spectral radii of the correlation matrix of the internal units and of the conceptor matrices
settings.plot_err_distr = true; % True = plot the distribution of the NRMSE over the specified number of runs
settings.investigate_pca = true; % True = investigate the set_precision_pca alg, comapring it to set_precision_rand or set_precision_randunif
settings.investigate_pca_rand = 'set_precision_rand'; % possible options: 'set_precision_rand' or 'set_precision_randunif'
settings.verbose = true; % for exhaustive search
settings.no_runs = 3;

%% Run algorithm
if settings.investigate_pca
    settings.no_runs = 3;
    maxdvec = [];
    maxdvecRand = [];
    uvvecPCA = [];
    uvvecRand = [];
    NRMSEvec = zeros(1, settings.no_runs);
    NRMSEvec_rand = zeros(1, settings.no_runs);
    settings.maxdvec = maxdvec;   
   
    settings.red_prec_alg = 'set_precision_pca';
    for i = 1:settings.no_runs
        tic;
        fprintf('run #%g/%g\n', i, settings.no_runs);
        output =  main_matrixC(settings);
        NRMSEvec(i) = output.meanNRMSE;
        maxdvec = [maxdvec, output.maxd];
        uvvecPCA = [uvvecPCA, output.uv];
        
        elapsed_time = toc;
        if i == 1
            est_run_time = elapsed_time*settings.no_runs*2;
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
    
    if strcmp(settings.investigate_pca_rand, 'set_precision_rand')
        settings.red_prec_alg = 'set_precision_rand';
    else
        settings.red_prec_alg = 'set_precision_randunif';
    end
    settings.maxdvec = maxdvec;
    for i = 1:settings.no_runs
        tic;
        fprintf('run #%g/%g\n', i, settings.no_runs);
        output =  main_matrixC(settings);
        NRMSEvec_rand(i) = output.meanNRMSE;
        maxdvecRand = [maxdvecRand, output.maxd];
        uvvecRand = [uvvecRand, output.uv];
        
        elapsed_time = toc;
        est_run_time = est_run_time - elapsed_time;
        if est_run_time < 60
            fprintf('Estimated remaining run-time: %g seconds \n', est_run_time);
        else
            fprintf('Estimated remaining run-time: %g minutes \n', est_run_time/60);
        end
    end
    plot_pca_investigations(maxdvec, maxdvecRand, uvvecPCA, uvvecRand, NRMSEvec, NRMSEvec_rand, settings.no_runs);
else
    % if we don't investigate the pca algorithm
    NRMSEvec = zeros(1, settings.no_runs);
    for i = 1:settings.no_runs
        tic;
        fprintf('run #%g/%g\n', i, settings.no_runs);
        output =  main_matrixC(settings);
        NRMSEvec(i) = output.meanNRMSE;
        
        elapsed_time = toc;
        if i == 1
            est_run_time = elapsed_time*settings.no_runs;
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
    
    fprintf('Results for %g runs, precision %g:\n', settings.no_runs, 2*settings.b + 1);
    fprintf('Mean meanNRMSE %g\n', mean(NRMSEvec));
    fprintf('Min meanNRMSE %g\n', min(NRMSEvec));
    fprintf('Max meanNRMSE %g\n', max(NRMSEvec));
    fprintf('Variance meanNRMSE %g\n', var(NRMSEvec));
    fprintf('Standard deviation meanNRMSE %g\n', std(NRMSEvec));
    
    % look at the distribution of testing errors
    if settings.plot_err_distr
        f = figure(); clf;
        set(gca,'fontsize',20)
        hold on;
        hist(NRMSEvec, 30)
        title(sprintf('Distribution NRMSE test\n mean = %g, var = %g, std = %g', mean(NRMSEvec), var(NRMSEvec), std(NRMSEvec)))
        xlabel('values')
        ylabel('counts')
        saveas(f, 'images/NRMSEtest.png', 'png')
        saveas(f, 'images/NRMSEtest.fig', 'fig')
        hold off;
    end
end
end