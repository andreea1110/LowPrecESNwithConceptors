function plot_patterns(patternsInt, netOutputInt, t, plotLength, meanNRMSE, train_xPL, neuronNr1, neuronNr2, SRcollector, C, settings)
%ploting the reference signal versus the ESN output
k = size(patternsInt, 1);
if settings.plot_full
    k = size(patternsInt, 1);
    for i = 1:k
        figure();
        clf;
        subplot(2, 2, 1);
        plot(t, patternsInt(i, 1:size(t, 2)), 'b', 'Linewidth', 2);
        hold on;
        plot(t, netOutputInt(i, 1:size(t, 2)), 'r--', 'Linewidth', 2);
        xlim([1, plotLength]);
        legend('reference', 'ESN output');
        title('Reference pattern vs. ESN output');
        
        subplot(2, 2, 2);
        plot(t, train_xPL{1, i}(neuronNr1, :), 'g', 'Linewidth', 2);
        hold on;
        plot(t, train_xPL{1, i}(neuronNr2, :), 'k', 'Linewidth', 2);
        hold off;
        title('Two neurons')
        
        %plotting the log_10  of the singular values of R, the correlation
        %matrix of the internal unit's activations
        subplot(2, 2, 3);
        hold on;
        plot(log10(diag(SRcollector{i})),'k','LineWidth',3);
        plot(zeros(1,length(diag(SRcollector{i}))),'k--');
        hold off;
        title('log_{10}(\sigma)');
        
        %plotting the log_10  of the singular values of C, the conceptor
        %matrix
        subplot(2, 2, 4);
        hold on;
        plot(diag(C.s{i}),'k','LineWidth',3);
        %plot(zeros(1,length(C.s{i})),'k--');
        %plot(ones(1,length(C.s{i})),'k--');
        hold off;
        title('s (singular values of C)')
        
    end
else
    f = figure(); clf;
    for i = 1:k
        subplot(2, 2, i);
        set(gca,'fontsize',10);
        hold on;
        plot(t, patternsInt(i, 1:size(t, 2)), 'b', 'Linewidth', 2);
        hold on;
        plot(t, netOutputInt(i, 1:size(t, 2)), 'r--', 'Linewidth', 2);
        %title(sprintf('NRMSE: %g', NRMSE_vec(i)));
        xlim([1, plotLength]);
        xlabel('Time [time units]')
        ylabel('Signal [signal units]')
        legend('reference', 'ESN output');
        supertitle(sprintf('Mean NRMSE: %g', meanNRMSE));
    end
    saveas(f, 'images/signals.png', 'png')
    saveas(f, 'images/signals.fig', 'fig')
end
end