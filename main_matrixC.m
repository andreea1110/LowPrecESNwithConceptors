function output = main_matrixC(settings)
%% Description
%%% Training a ESN with conceptors to store dynamical patterns, using the
%%% original matrix-conceptor algorithm.
%%% The network has 1 input neuron and 1 output neuron.
%%% The number of neurons in the dynamical reservoir is adjustable.

%% Setting ESN parameters
ESN.DR.size = 100; % size of the dynamical reservoir (DR)
ESN.DR.SR = 1.5; % size of the spectral radius (SR) of the weight matrix W of the DR
ESN.DR.connectivity = 1; % connectivity of the units in the DR -> reccomended: 1 for reduced precision and 0.1 for full precision
ESN.inputScaling = 1.5; % scaling of the input-to-DR weights (pattern feeding weights) %1.5 initially
ESN.biasScaling = 0.2; % scaling of the bias vector
ESN.regW  = 0.0001; % regularizer used in the ridge (Tychonov) regression, for training W
ESN.washoutLength = 500; % number of discrete time steps disregarded before training; the time required for the ESN to wash-out the initial transient dynamics
ESN.learnLength = 1000; % size of the training set, here discrete signal points
ESN.testLength = 200; % size of the testing set
ESN.plotLength = 20; % length of the plotted signal
ESN.regWout  = 0.01; % regularizer used in the ridge (Tychonov) regression, for training the DR-to-output weights W^{out}
ESN.WparameterNL = 0.0; % ?
ESN.stateNL = 0.0; % ?
ESN.noiselevel = 0.000;
ESN.WparameterNL = 0.0; % ?
intRate = 20; % interpolation rate for the signals

%% Setting Conceptor parameters
%C.reg = 0.6; % regularizer used in the ridge (Tikhonov) regression, for training the conceptor matrix
C.svPlotLength = 50; % singular value (of the conceptor matrix) plot length?
C.aperture = 10; % used as a regularizer for the regression problem solved when computing the conceptor matrices

%% Select patterns from list
patterns = [53 54 10 36];

%% Train: store patterns into ESN with Conceptors
[ESN, C, NRMSE_W, NRMSE_readout, train_xPL, train_pPL, SRcollector, maxd1, uv1] = trainESNC(ESN, C, patterns, settings);
%add_noise_to_W(ESN); % ?

%% Reduce precision of all trained parameters

if settings.red_prec % if we have to reduce the precision
    % select the function which will be used in order to reduce the precision
    % of the parameters, according to the specified settings
    if strcmp(settings.red_prec_alg, 'set_precision_orig')
        set_precision = @set_precision_orig;
    elseif strcmp(settings.red_prec_alg, 'set_precision_distr1')
        set_precision = @set_precision_distr1;
    elseif strcmp(settings.red_prec_alg, 'set_precision_distr2')
        set_precision = @set_precision_distr2;
    elseif strcmp(settings.red_prec_alg, 'set_precision_pca')
        set_precision = @set_precision_pca;
    elseif strcmp(settings.red_prec_alg, 'set_precision_sc')
        set_precision = @set_precision_sc;
    elseif strcmp(settings.red_prec_alg, 'set_precision_rand')
        set_precision = @set_precision_rand;
        parameters.maxdvec = settings.maxdvec;
    else
        set_precision = @set_precision_exsearch;
    end
    
    if settings.plot_init_distr
        % plotting the distribution of the values in the Wout and W matrices
        figure();
        hist(ESN.Wout(:), 40);
        title('ESN.Wout');
        xlabel('values')
        ylabel('counts')
        
        figure();
        hist(ESN.W(:), 40);
        title('ESN.W');
        xlabel('values')
        ylabel('counts')
        
        figure() % used for plotting the distribution of the C matrices
    end
    
    if settings.investigate_pca
        maxd2 = zeros(1, 2 + length(patterns));
        uv2 = zeros(1, 2 + length(patterns));
        parameters.maxdvec = settings.maxdvec;
    end
    
    % set general parameters for reducing the precision
    parameters.ESN = ESN;
    parameters.C = C;
    parameters.patterns = patterns;
    parameters.train_pPL = train_pPL;
    parameters.intRate = intRate;
    parameters.b = settings.b;
    
    % reduce the precision of Wout
    parameters.D = ESN.Wout;
    parameters.mat = 'Wout';
    res = set_precision(parameters);
    ESN.Wout = res.D_lp;
    if settings.investigate_pca
        maxd2(1, 1) = res.maxd;
        uv2(1, 1) = res.uv;
    end
    
    % reduce the precision of W
    parameters.D = ESN.W;
    parameters.mat = 'W';
    res = set_precision(parameters);
    ESN.W = res.D_lp;
    if settings.investigate_pca
        maxd2(1, 2) = res.maxd;
        uv2(1, 2) = res.uv;
    end
    
    
    for i = 1:length(patterns)
        if settings.plot_init_distr
            subplot(220 + i)
            hist(C.mat{i}(:), 40);
            title(sprintf('Conceptor matrix %g', i));
            xlabel('values');
            ylabel('counts');
        end
        % reduce the precision of Ci
        parameters.D = C.mat{i};
        parameters.mat = strcat('C', int2str(i));
        res = set_precision(parameters);
        C.mat{i} = res.D_lp;
        if settings.investigate_pca
            maxd2(1, 2 + i) = res.maxd;
            uv2(1, 2 + i) = res.uv;
        end
    end
    
    if settings.investigate_pca
        output.maxd = [maxd1, maxd2];
        output.uv = [uv1, uv2];
    end
end
%% Test
netOutput = testESNC(ESN, C, patterns);
[netOutputInt, patternsInt] = interpolate(netOutput, train_pPL, intRate, ESN);
[netOutputIntAl, patternsIntResized] = align_signals(netOutputInt, patternsInt, intRate, ESN);
%fprintf('mean NRMSE W: %0.3g\n', mean(NRMSE_W));
%fprintf('NRMSE readout (training error): %g\n', NRMSE_readout);
output.meanNRMSE = compute_error(netOutputIntAl, patternsIntResized, size(patternsIntResized, 2));

%% Plotting
if settings.plot_signals
    t = 1:ESN.plotLength;
    neuronNr1 = 71; neuronNr2 = 80;
    plot_patterns(patternsIntResized, netOutputIntAl, t, ESN.plotLength, output.meanNRMSE, train_xPL, neuronNr1, neuronNr2, SRcollector, C, settings);
end
end
