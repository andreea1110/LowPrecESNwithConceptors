%%%% plain demo that when a RNN is driven by different signals, the induced
%%%% internal signals will inhabit different subspaces of the signal space.
%%%% here we run this basic conceptor demo by first training it as usual,
%%%% then reducing the bit accuracy of all weights

set(0,'DefaultFigureWindowStyle','docked');

%%% Experiment control
randstate = 1; newNets = 1; newSystemScalings = 1;
pattHandles;

%%% Setting system params
Netsize = 100; % network size
NetSR = 1.5; % spectral radius
NetinpScaling = 1.5; % scaling of pattern feeding weights
BiasScaling = 0.2; % size of bias
precFixedWeights = 2; interpFixed = 0; % precision of all parameters. Set to Inf if no truncation is done.
precWout = Inf; interpWout = 0;
precW = 2; interpW = 1;
precC = 5; interpC = 1;

%%% loading learning
TychonovAlpha = .0001; % regularizer for  W training
washoutLength = 500;
learnLength = 1000;
signalPlotLength = 20;

%%% pattern readout learning
TychonovAlphaReadout = .01;

%%% C learning and testing
alpha = 10;
CtestLength = 200;
SplotLength = 50;

WparameterNL = 0.0;
stateNL = 0.0;


%%% Setting patterns

patterns = [53 54 10 36];
%patterns = [23 6];
%patterns = [1 2 21 20 22 8 19 6  16 9 10 11 12];

% 1 = sine10  2 = sine15  3 = sine20  4 = spike20
% 5 = spike10 6 = spike7  7 = 0   8 = 1
% 9 = rand5; 10 = rand5  11 = rand6 12 = rand7
% 13 = rand8 14 = sine10range01 15 = sine10rangept5pt9
% 16 = rand3 17 = rand9 18 = rand10 19 = 0.8 20 = sineroot27
% 21 = sineroot19 22 = sineroot50 23 = sineroot75
% 24 = sineroot10 25 = sineroot110 26 = sineroot75tenth
% 27 = sineroots20plus40  28 = sineroot75third
% 29 = sineroot243  30 = sineroot150  31 = sineroot200
% 32 = sine10.587352723 33 = sine10.387352723
% 34 = rand7  35 = sine12  36 = 10+perturb  37 = sine11
% 38 = sine10.17352723  39 = sine5 40 = sine6
% 41 = sine7 42 = sine8  43 = sine9 44 = sine12
% 45 = sine13  46 = sine14  47 = sine10.8342522
% 48 = sine11.8342522  49 = sine12.8342522  50 = sine13.1900453
% 51 = sine7.1900453  52 = sine7.8342522  53 = sine8.8342522
% 54 = sine9.8342522 55 = sine5.19004  56 = sine5.8045
% 57 = sine6.49004 58 = sine6.9004 59 = sine13.9004
% 60 = 18+perturb

%%% Initializations

%randn('state', randstate);
%rand('twister', randstate);

% Create raw weights
if newNets
    if Netsize <= 20
        Netconnectivity = 1;
    else
        Netconnectivity = 1; %10/Netsize;
    end
    WstarRaw = generate_internal_weights(Netsize, Netconnectivity);
    WinRaw = rand(Netsize, 1) - 0.5;
    WbiasRaw = rand(Netsize, 1) - 0.5;
end

% Scale raw weights and initialize weights
if newSystemScalings
    Wstar = NetSR * WstarRaw;
    Win = NetinpScaling * WinRaw;
    Wbias = BiasScaling * WbiasRaw;
end

% Set pattern handles
%pattHandles; 
%%% in order to get same patterns as in the hierarchical
%architecture demo, run that demo first and do not call pattHandles here
%again. 

I = eye(Netsize);

% % learn equi weights

% harvest data from network externally driven by patterns
Np = length(patterns);
allTrainArgs = zeros(Netsize, Np * learnLength);
allTrainOldArgs = zeros(Netsize, Np * learnLength);
%allTrainTargs = zeros(Netsize, Np * learnLength);
allTrainOuts = zeros(1, Np * learnLength);
readoutWeights = cell(1,Np);
%patternCollectors = cell(1,Np);
%xCollectorsCentered = cell(1,Np);
%xCollectors = cell(1,Np);
SRCollectors = cell(1,Np);
%URCollectors = cell(1,Np);
patternRs = cell(1,Np);
train_xPL = cell(1,Np);
train_pPL = cell(1,Np);
%startXs = zeros(Netsize, Np);
% collect data from driving native reservoir with different drivers
for p = 1:Np % for each of the input patterns
    patt = patts{patterns(p)}; % current pattern generator
    xCollector = zeros(Netsize, learnLength );
    xOldCollector = zeros(Netsize, learnLength );
    pCollector = zeros(1, learnLength );
    x = zeros(Netsize, 1);
    for n = 1:(washoutLength + learnLength)
        u = patt(n); % pattern input
        xOld = x;
        x = tanh(Wstar * x + Win * u + Wbias);
        if n > washoutLength
            xCollector(:, n - washoutLength ) = x;
            xOldCollector(:, n - washoutLength ) = xOld;
            pCollector(1, n - washoutLength) = u;
        end
    end
    
    %xCollectorCentered = xCollector - ...
    %    repmat( mean(xCollector,2),1,learnLength);
    %xCollectorsCentered{1,p} = xCollectorCentered;
    %xCollectors{1,p} = xCollector;
    R = xCollector * xCollector' / learnLength;
    [Ux Sx Vx] = svd(R);
    SRCollectors{1, p} = Sx;
    %URCollectors{1,p} = Ux;
    patternRs{p} = R;
    
    
    %startXs(:,p) = x;
    train_xPL{1,p} = xCollector(:,1:signalPlotLength);
    train_pPL{1,p} = pCollector(1,1:signalPlotLength);
    
    %patternCollectors{1,p} = pCollector;
    allTrainArgs(:, (p-1)*learnLength+1:p*learnLength) = ...
        xCollector;
    allTrainOldArgs(:, (p-1)*learnLength+1:p*learnLength) = ...
        xOldCollector;
    allTrainOuts(1, (p-1)*learnLength+1:p*learnLength) = ...
        pCollector;
    %allTrainTargs(:, (p-1)*learnLength+1:p*learnLength) = ...
        %Win * pCollector;
end

%%% compute readout

Wout = ((allTrainArgs * allTrainArgs' + ...
    TychonovAlphaReadout * I)\allTrainArgs * allTrainOuts')';
% training error
NRMSE_readout = nrmse(Wout*allTrainArgs, allTrainOuts);


%%% compute W
Wtargets = (atanh(allTrainArgs) - repmat(Wbias,1,Np*learnLength));
allTrainOldArgs2 = [zeros(Netsize, 1), allTrainArgs];
allTrainOldArgs2(:, end) = [];
W = ((allTrainOldArgs2 * allTrainOldArgs2' + ...
    TychonovAlpha * I)\allTrainOldArgs2 * Wtargets')';
% training errors per neuron
NRMSE_W = nrmse(W*allTrainOldArgs, Wtargets);

% % compute conceptors
Cs = cell(4, Np);
aperture = 10;
for p = 1:Np
    R = patternRs{p};
    
    Cs{1, p} = R / (R + (aperture)^(-2) * I);
    
    [U S V] = svd(R);
    Snew = (S /(S + alpha^(-2) * eye(Netsize)));
    
    C= R / (R + (aperture)^(-2) * I);
    
    %C = U * Snew * U';
    Cs{1, p} = C;
    %Cs{2, p} = U;
    Cs{3, p} = diag(Snew);
    %Cs{4, p} = diag(S);
end

% % test with C
%%
% add noise to W
Wnoisy = W + abs(W) .* (WparameterNL * (rand(Netsize,Netsize)-0.5));

%x_CTestPL = zeros(5, CtestLength, Np);
p_CTestPL = zeros(1, CtestLength, Np);
for p = 1:Np
    C = Cs{1, p};
    %x = startXs(:,p);
    x = 0.5*randn(Netsize,1);
    
    for n = 1:CtestLength
        x = tanh(Wnoisy *  x + Wbias )+ stateNL * (rand(Netsize,1)-0.5);
        x = C * x;
        %x_CTestPL(:,n,p) = x(1:5,1);
        p_CTestPL(:,n,p) = Wout * x;
    end
end



%% plotting


test_pAligned_PL = cell(1,Np);
%test_xAligned_PL = cell(1,Np);
NRMSEsAligned = zeros(1,Np);
%MSEsAligned = zeros(1,Np);

for p = 1:Np
    intRate = 20;
    thisDriver = train_pPL{1,p};
    thisOut = p_CTestPL(1,:,p);
    thisDriverInt = interp1((1:signalPlotLength)',thisDriver', (1:(1/intRate):signalPlotLength)', 'spline')';
    thisOutInt = interp1((1:CtestLength)', thisOut', (1:(1/intRate):CtestLength)', 'spline')';
    
    L = size(thisOutInt,2); 
    M = size(thisDriverInt,2);
    phasematches = zeros(1,L - M);
    for phaseshift = 1:(L - M)
        phasematches(1,phaseshift) = norm(thisDriverInt - thisOutInt(phaseshift:phaseshift+M-1));
    end
    %[maxVal maxInd] = max(-phasematches);
    [maxVal maxInd] = min(phasematches);
    test_pAligned_PL{1,p} = thisOutInt(maxInd:intRate:(maxInd+intRate*signalPlotLength-1));
    %coarseMaxInd = ceil(maxInd / intRate);
    %test_xAligned_PL{1,p} = x_CTestPL(:,coarseMaxInd:coarseMaxInd+signalPlotLength-1, p);
    NRMSEsAligned(1,p) = nrmse(test_pAligned_PL{1,p},train_pPL{1,p});
    %MSEsAligned(1,p) = mean((test_pAligned_PL{1,p} - train_pPL{1,p}).^2);
end
disp('******************');
fprintf('NRMSE readout: %g\n', NRMSE_readout);
fprintf('mean NRMSE W: %0.3g\n  meanAbsW: %0.3g\n', mean(NRMSE_W), mean(mean(abs(W))));
disp(['NRMSEsAligned = ' num2str(NRMSEsAligned, ' %0.2g\t')]);
fprintf('meanNRMSE = %0.3g\n', mean(NRMSEsAligned));

%%
figure(2); clf;
fs = 10; fsNRMSE = 10;
%set(gcf,'DefaultAxesColorOrder',[0  0.4 0.65 0.8]'*[1 1 1]);
set(gcf, 'WindowStyle','normal');

set(gcf,'Position', [850 400 800 500]);

pick1 = 71; pick2 = 80;
col1 = 0.6*[1 1 1]; col2 = 0.3*[1 1 1];
for p = 1:Np
    subplot(Np,4,(p-1)*4+1);
    plot(test_pAligned_PL{1,p}, 'LineWidth',10, 'Color',0.75*[1 1 1]); hold on;
    plot(train_pPL{1,p},'k','LineWidth',1.5); hold off;
    if p == 1
        title('p and y','FontSize',fs);
    end
    if p ~= Np
        set(gca, 'XTickLabel',[]);
    end
    set(gca, 'YLim',[-1,1], 'FontSize',fs);
%     rectangle('Position', [0.5,-0.95,8,0.5],'FaceColor','w',...
%         'LineWidth',1);
%      text(1,-0.7,num2str(NRMSEsAligned(1,p),2),...
%         'Color','k','FontSize',fsNRMSE, 'FontWeight', 'bold');
    
    subplot(Np,4,(p-1)*4+2); hold on;
    plot(train_xPL{1,p}(pick1,:)','Color',col1,'LineWidth',3);
    plot(train_xPL{1,p}(pick2,:)','Color',col2,'LineWidth',3);
    hold off;
    if p == 1
        title('two neurons','FontSize',fs);
    end
    if p ~= Np
        set(gca, 'XTickLabel',[]);
    end
    set(gca,'YLim',[-1,1], 'FontSize',fs, 'Box', 'on');
    
    subplot(Np,4,(p-1)*4+3);
    %diagNormalized = sDiagCollectors{1,p} / sum(sDiagCollectors{1,p});
    hold on;
    plot(log10(diag(SRCollectors{1,p})),'k','LineWidth',3);
    plot(zeros(1,Netsize),'k--');
    hold off;
    
    set(gca,'YLim',[-17,5], 'YTick',[-10, 0], 'FontSize',fs, 'Box','on');
    if p == 1
        title('log10 \sigma','FontSize',fs);
    end
    if p < 4
        set(gca,'XTick',[]);
    end
    
    subplot(Np,4,(p-1)*4+4);
    hold on;
    plot(zeros(1,Netsize),'k--');
    plot(ones(1,Netsize),'k--');
    plot(Cs{3, p},'k','LineWidth',3);
    hold off;
    if p == 1
        title('s','FontSize',fs);
    end
    set(gca,'YLim',[-0.1,1.1], 'YTick',[0 1], 'FontSize',fs, 'Box','on');
    if p < 4
        set(gca,'XTick',[]);
    end
end