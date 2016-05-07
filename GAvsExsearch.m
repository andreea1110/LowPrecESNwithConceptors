clear all; 
clc;

nrmse_exsearch = load('variables\GAvsExsearch\NRMSE_set_precision_exsearch.mat');
nrmse_ga = load('variables\GAvsExsearch\NRMSE_set_precision_ga.mat');

f = figure(); clf;
set(gca,'fontsize',16);
hold on;
%map = brewermap(3,'Set1'); 
histf(nrmse_exsearch.NRMSEvec, 0:.1:2, 'facecolor', 'b','facealpha',0.5, 'edgecolor','none');
hold on;
histf(nrmse_ga.NRMSEvec, 0:.1:2, 'facecolor', 'r','facealpha',0.5, 'edgecolor','none');
xlabel('values')
ylabel('counts')
title(sprintf('Error (NRMSE) distrubution for 10 runs of \n the genetic algorithm and the exhaustive search algorithm'));
legalpha('exhaustive search','genetic algorithm')
%savefig('results/histograms/PCAvsRand.fig')
saveas(f, 'images/GAvsExsearch.png', 'png')
saveas(f, 'images/GAvsExsearch.fig', 'fig')