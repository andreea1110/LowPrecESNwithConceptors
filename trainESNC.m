function [ESN, C, NRMSE_W, NRMSE_readout, train_xPL, train_pPL, SRcollector, maxd1, uv1] = trainESNC(ESN, C, patterns, settings)
% Stores the signals, stored as function pointers in the vector patterns,
% in the DR of the ESN.
%
% Notations:
%   W = matrix storing the weights of the connections among the internal
%      neurons, after training
%   Win = matrix storing the weights of the connections among the input and
%       internal neurons
%   Wout = matrix storing the weights of the connections among the internal
%        and output neurons, after training
%
%   u = value of the input unit
%   x = activations of the internal units
%   y = value of the ouput unit
%

%% Initializations
k = length(patterns); % number of patterns to be stored in the ESN
x_AllTrain = zeros(ESN.DR.size, k*ESN.learnLength); % allTrainArgs
y_AllTrain = zeros(1, k*ESN.learnLength); % allTrainOuts
train_xPL = cell(1, k);
train_pPL = cell(1, k);
Rcollector = cell(1, k); % collects the correlation matrices for each input pattern
SRcollector = cell(1, k); % collects the singular value matrices S, for the correlation matrix of the internal states R, corresponding to each pattern
I = eye(ESN.DR.size); % identity matrix of the size of W

% Randomly initializing the weights of the dynamical reservoir neuron
% connections, the input weights and the weights of the bias vector
WstarRaw = generate_internal_weights(ESN.DR.size, ESN.DR.connectivity);
WinRaw = rand(ESN.DR.size, 1) - 0.5;
WbiasRaw = rand(ESN.DR.size, 1) - 0.5;


% Scaling the previously initialized weights
ESN.Wstar = ESN.DR.SR * WstarRaw;
ESN.Win = ESN.inputScaling * WinRaw;
ESN.Wbias = ESN.biasScaling * WbiasRaw;

%% Reduce the precision of the randomly initialized weight matrices
maxd1 = zeros(1, 3);
uv1 = zeros(1, 3);

if settings.red_prec % if we have to reduce the precision
    % select the function which will be used in order to reduce the precision
    % of the parameters, according to the specified settings
    if strcmp(settings.red_prec_alg, 'set_precision_orig')
        set_precision = @set_precision_orig;
    elseif strcmp(settings.red_prec_alg, 'set_precision_distr1')
        set_precision = @set_precision_distr1;
    elseif strcmp(settings.red_prec_alg, 'set_precision_pca')
        set_precision = @set_precision_pca;
    elseif strcmp(settings.red_prec_alg, 'set_precision_rand')
        set_precision = @set_precision_rand;
    elseif strcmp(settings.red_prec_alg, 'set_precision_randunif')
        set_precision = @set_precision_randunif;
    elseif strcmp(settings.red_prec_alg, 'set_precision_sc')
        set_precision = @set_precision_sc;
    else
        set_precision = @set_precision_id; % identity funtion => don't reduce the precision of those matrices
    end
    
    if settings.investigate_pca
        parameters.maxdvec = settings.maxdvec;
    end
    
    % reducing the precision of the matrices
    parameters.b = settings.b;
    
    settings.svd = false;
    if settings.svd
        % reduce the precision of the svd of Win
        [UWin SWin VWin] = svd(ESN.Win);
        parameters.D = UWin;
        res = set_precision(parameters);
        UWin = res.D_lp;
        parameters.D = SWin;
        res = set_precision(parameters);
        SWin = res.D_lp;
        parameters.D = VWin;
        res = set_precision(parameters);
        VWin = res.D_lp;
        ESN.Win = UWin*SWin*VWin';
        
        % reduce the precision of the svd of Wbias
        [UWbias SWbias VWbias] = svd(ESN.Wbias);
        parameters.D = UWbias;
        res = set_precision(parameters);
        UWbias = res.D_lp;
        parameters.D = SWbias;
        res = set_precision(parameters);
        SWbias = res.D_lp;
        parameters.D = VWbias;
        res = set_precision(parameters);
        VWbias = res.D_lp;
        ESN.Wbias = UWbias*SWbias*VWbias';
        
        % reduce the precision of the svd of Wstar
        [UWstar SWstar VWstar] = svds(ESN.Wstar);
        parameters.D = UWstar;
        res = set_precision(parameters);
        UWstar = res.D_lp;
        parameters.D = SWstar;
        res = set_precision(parameters);
        SWstar = res.D_lp;
        parameters.D = VWstar;
        res = set_precision(parameters);
        VWstar = res.D_lp;
        ESN.Wstar = UWstar*SWstar*VWstar';
        
    else
        % reduce the precision of Win
        parameters.D = ESN.Win;
        res = set_precision(parameters);
        ESN.Win = res.D_lp;
        if settings.investigate_pca
            maxd1(1, 1) = res.maxd;
            uv1(1, 1) = res.uv;
        end
        
        % reduce the precision of Wbias
        parameters.D = ESN.Wbias;
        res = set_precision(parameters);
        ESN.Wbias = res.D_lp;
        if settings.investigate_pca
            maxd1(1, 2) = res.maxd;
            uv1(1, 2) = res.uv;
        end
        
        % reduce the precision of Wstar
        parameters.D = ESN.Wstar;
        res = set_precision(parameters);
        ESN.Wstar = res.D_lp;
        if settings.investigate_pca
            maxd1(1, 3) = res.maxd;
            uv1(1, 3) = res.uv;
        end
    end
end
%% Harvest data from network externally driven by the patterns
pattHandles;
% initialize noise vector
noise = rand(1, ESN.washoutLength + ESN.learnLength) - 0.5;


for i = 1:k % for each of the input patterns
    u = patts{patterns(i)}; % current pattern generator
    xCollector = zeros(ESN.DR.size, ESN.learnLength); %xCollector
    pCollector = zeros(1, ESN.learnLength); % pattern collector
    x = zeros(ESN.DR.size, 1); % initialize the activations of the internal units to 0
    for j = 1:(ESN.washoutLength + ESN.learnLength)
        x = tanh(ESN.Wstar * x + ESN.Win * u(j) + ESN.Wbias) + ESN.noiselevel * 2.0 * noise(j); % compute the activations of the internal units at time step j
        if j > ESN.washoutLength
            xCollector(:, j - ESN.washoutLength) = x;
            pCollector(1, j - ESN.washoutLength) = u(j);
        end
    end
    
    % computing the correlation matrix of the internal units' activations,
    % for pattern i
    R = (xCollector * xCollector')/ESN.learnLength;
    Rcollector{i} = R;
    
    % computing the singular value decomposition of R; Ux*Sx*Vx' = R;
    [Ux Sx Vx] = svd(R);
    SRcollector{i} = Sx;
    
    train_xPL{1,i} = xCollector(:,1:ESN.plotLength); % why?
    train_pPL{1,i} = pCollector(1,1:ESN.plotLength); % why?
    
    x_AllTrain(:, (i-1)*ESN.learnLength+1:i*ESN.learnLength) = xCollector;
    y_AllTrain(1, (i-1)*ESN.learnLength+1:i*ESN.learnLength) = pCollector;
end

%% Compute output weights using ridge regression
ESN.Wout = ((x_AllTrain * x_AllTrain' + ESN.regWout * I)\x_AllTrain * y_AllTrain')';
% training error
NRMSE_readout = nrmse(ESN.Wout*x_AllTrain, y_AllTrain);

%% Recomputing W (a.k.a. loading the DR with patterns)
Wtargets = atanh(x_AllTrain) - repmat(ESN.Wbias,1,k*ESN.learnLength);
x_AllTrainOld = [zeros(ESN.DR.size, 1), x_AllTrain]; % ???
x_AllTrainOld(:, end) = [];
ESN.W = ((x_AllTrainOld * x_AllTrainOld' + ESN.regW * I)\x_AllTrainOld * Wtargets')';
% training errors per neuron
NRMSE_W = nrmse(ESN.W*x_AllTrainOld, Wtargets);

%% Computing the conceptor matrices
C.mat = cell(1, k);
C.s = cell(1, k);
for i = 1:k % iterating through the patterns
    R = Rcollector{i}; % retrieve the corresponding correlation matrix
    C.mat{i} = R / (R + (C.aperture)^(-2) * I);
    [U S V] = svd(C.mat{i});
    C.s{i} = S;
end
end
