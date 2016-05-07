function plot_pca_investigations(maxdvec, maxdvecRand, uvvecPCA, uvvecRand, NRMSEvec, NRMSEvec_rand, noRuns)
fontsize = 16;
% look at the error compared to number of values in the matrices
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
x = 1:length(NRMSEvec);
[AX,H1,H2] = plotyy(x, NRMSEvec, x, sum(reshape(uvvecPCA, 9, noRuns), 1),'plot');
%disp(cov(NRMSEvec, sum(reshape(uvvecPCA, 9, noRuns), 1)));
c = cov(NRMSEvec, sum(reshape(uvvecPCA, 9, noRuns), 1));
set(AX,'fontsize',fontsize);
set(get(AX(1),'Ylabel'),'String','NRMSE','fontsize',fontsize) 
set(get(AX(2),'Ylabel'),'String','sum no. uniq. el', 'fontsize',fontsize, 'Color', 'r')
set(AX(2), 'YColor', 'r');
set(H1,'Color', 'b', 'Linewidth', 2)
set(H2,'Color', 'r', 'LineStyle','--','Linewidth', 2)
title(strcat('NRMSE compared to number of values in the matrices,', sprintf(' cov = %g.', c(2, 1))));
xlabel('runs')
%legend('NRMSE', 'sum no. uniq. el')
saveas(f, 'images/errNoVals.png', 'png')
saveas(f, 'images/errNoVals.fig', 'fig')

% plot error vs. average precision from each matrix precision 
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
scatter(NRMSEvec, log2(sum(reshape(uvvecPCA, 9, noRuns), 1)/9) ,'filled');
xlabel('NRMSE');
ylabel('average precision [bits]');
c = cov(NRMSEvec, log2(sum(reshape(uvvecPCA, 9, noRuns), 1)/9));
title(sprintf('Precision vs. NRMSE for %g runs, cov = %g', noRuns, c(2, 1)));
saveas(f, 'images/errorVsPrecision.png', 'png');
saveas(f, 'images/errorVsPrecision.fig', 'fig');


% look at the distribution of the maxd values for pca
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
hist(maxdvec, 30)
title(sprintf('maxd distribution PCA\n mean = %g, var = %g, std = %g', mean(maxdvec), var(maxdvec), std(maxdvec)));
xlabel('values')
ylabel('counts')
%savefig('results/histograms/maxd.fig')
saveas(f, 'images/maxdPCA.png', 'png')
saveas(f, 'images/maxdPCA.fig', 'fig')

% look at the distribution of the maxd values for rand
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
hist(maxdvecRand, 30)
title(sprintf('maxd distribution rand\n mean = %g, var = %g, std = %g', mean(maxdvecRand), var(maxdvecRand), std(maxdvecRand)))
xlabel('values')
ylabel('counts')
%savefig('results/histograms/maxd.fig')
saveas(f, 'images/maxdRand.png', 'png')
saveas(f, 'images/maxdRand.fig', 'fig')

% look at the distribution of the number of unique values in the matrices
% for the PCA algorithm
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
hist(uvvecPCA, 30)
title(sprintf('number of unique values in the reduced-precision matrices (PCA)\n mean = %g, var = %g, std = %g', mean(uvvecPCA), var(uvvecPCA), std(uvvecPCA)))
xlabel('values')
ylabel('counts')
%savefig('results/histograms/uvmatPCA.fig')
saveas(f, 'images/uvmatPCA.png', 'png')
saveas(f, 'images/uvmatPCA.fig', 'fig')


% look at the distribution of the number of unique values in the matrices
% for the random algorithm
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
hist(uvvecRand, 30)
title(sprintf('number of unique values in the reduced-precision matrices (rand)\n mean = %g, var = %g, std = %g', mean(uvvecRand), var(uvvecRand), std(uvvecRand)))
xlabel('values')
ylabel('counts')
%savefig('results/histograms/uvmatRand.fig')
saveas(f, 'images/uvmatRand.png', 'png')
saveas(f, 'images/uvmatRand.fig', 'fig')

% look at the distribution of the error for the PCA algorithm
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
hist(NRMSEvec, 30)
title(sprintf('error (NRMSE) distrubution for the PCA algorithm\n mean = %g, var = %g, std = %g', mean(NRMSEvec), var(NRMSEvec), std(NRMSEvec)))
xlabel('values')
ylabel('counts')
%savefig('results/histograms/NRMSEvec.fig')
saveas(f, 'images/NRMSEvec.png', 'png')
saveas(f, 'images/NRMSEvec.fig', 'fig')

% look at the distribution of the error for the random algorithm
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
hist(NRMSEvec_rand, 30)
title(sprintf('error (NRMSE) distrubution for the random algorithm\n mean = %g, var = %g, std = %g', mean(NRMSEvec_rand), var(NRMSEvec_rand), std(NRMSEvec_rand)))
xlabel('values')
ylabel('counts')
%savefig('results/histograms/NRMSEvec_rand.fig')
saveas(f, 'images/NRMSEvec_rand.png', 'png')
saveas(f, 'images/NRMSEvec_rand.fig', 'fig')

% look at the distribution of the error for the PCA algorithm vs. error in
% the random algorithm
f = figure(); clf;
set(gca,'fontsize',fontsize);
hold on;
%map = brewermap(3,'Set1'); 
histf(NRMSEvec, 0:.1:2, 'facecolor', 'b','facealpha',0.5, 'edgecolor','none')
hold on
histf(NRMSEvec_rand, 0:.1:2, 'facecolor', 'r','facealpha',0.5, 'edgecolor','none')
title('error (NRMSE) distrubution')
xlabel('values')
ylabel('counts')
legalpha('PCA','random')
%savefig('results/histograms/PCAvsRand.fig')
saveas(f, 'images/PCAvsRand.png', 'png')
saveas(f, 'images/PCAvsRand.fig', 'fig')
end