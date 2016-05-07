clear all;
clc;
load variables\exsearch100\paramsWout.mat

fontsize = 16;
% plot error vs. psi and phi
fig = figure();clf;
subplot(1, 2, 1)
surf(psi_range, phi_range, error);
colormap(jet);
xlabel('\psi', 'fontsize',fontsize);
ylabel('\phi', 'fontsize',fontsize);
zlabel('NRMSE', 'fontsize',fontsize);
title('Error vs. optimization parameters for Wout', 'fontsize',16);
hold on;
set(gca,'fontsize',fontsize);
grid on;
subplot(1, 2, 2)
% contour plot
contourf(psi_range, phi_range, error);
colormap(jet);
cbr = colorbar('location','southoutside');
set(cbr,'ZTick',0:0.1:2)
xlabel('\psi', 'fontsize',fontsize);
ylabel('\phi', 'fontsize',fontsize);
title('Error contour plot for Wout', 'fontsize',16);
hold on;
set(gca,'fontsize',fontsize);
grid on;
saveas(fig, strcat('images/surfContWout.png'), 'png');
saveas(fig, strcat('images/surfContWout.fig'), 'fig');


xvec = [];
yvec = [];
zvec = [];
x = psi_range;
y = phi_range;

for i = 1:length(psi_range)
    for j = 1:length(phi_range)
        xvec = [xvec, psi_range(i)];
        yvec = [yvec, phi_range(j)];
        zvec = [zvec, error(i, j)];
    end
end


p = polyfitn([xvec',yvec'],zvec,3);
polyn2sympoly(p);
fprintf('RMSE - polynomial of order 3: %g\n', p.RMSE);

f = zeros(length(psi_range), length(phi_range));
for i = 1:length(x)
    for j = 1:length(y)
        f(i, j) = -84.002*x(i)^3 + 159.7558*x(i)^2*y(j) + 40.8388*x(i)^2 - 172.8048*x(i)*y(j)^2 - 35.8391*x(i)*y(j) - 2.5267*x(i) - 383.5434*y(j)^3 + 179.0336*y(j)^2 - 17.0945*y(j) + 1.1181;

    end
end

[min_err, idx_error] = min(error(:));
fprintf('min_err = %g\n', min_err);
[min_f, idx_f] = min(f(:));
fprintf('min_f3 = %g\n', min_f);
[m, n] = ind2sub(size(f),idx_f);
err_aprox = error(m, n);
rel_err = (abs(min_err - err_aprox)/min_err)*100;
fprintf('Relative error (in percentage) of polyn. of order 3\n to the exhaustive search minimum: %g \n', rel_err );

% plot error vs. psi and phi
fig = figure();clf;
subplot(1, 2, 1)
surf(x, y, f);
colormap(jet);
xlabel('\psi', 'fontsize',fontsize);
ylabel('\phi', 'fontsize',fontsize);
zlabel('NRMSE', 'fontsize',fontsize);
title('Error - order 3 approximation', 'fontsize',16);
hold on;
set(gca,'fontsize',fontsize);
grid on;
subplot(1, 2, 2)
% contour plot
contourf(x, y, f);
colormap(jet);
cbr = colorbar('location','southoutside');
set(cbr,'ZTick',0:0.1:1.5);
xlabel('\psi', 'fontsize',fontsize);
ylabel('\phi', 'fontsize',fontsize);
title('Error approximation contour plot for Wout', 'fontsize',16);
hold on;
set(gca,'fontsize',fontsize);
grid on;
saveas(fig, strcat('images/Woutaprox3.png'), 'png');
saveas(fig, strcat('images/Woutaprox3.fig'), 'fig');


p = polyfitn([xvec',yvec'],zvec,4);
polyn2sympoly(p);
fprintf('RMSE - polynomial of order 4: %g\n', p.RMSE);

f = zeros(length(psi_range), length(phi_range));
for i = 1:length(x)
    for j = 1:length(y)
        f(i, j) = -442.1552*x(i)^4 - 285.2892*x(i)^3*y(j) + 298.2511*x(i)^3 - 923.6518*x(i)^2*y(j)^2 + 515.6597*x(i)^2*y(j) - 73.0306*x(i)^2 + 3722.3316*x(i)*y(j)^3 - 920.0436*x(i)*y(j)^2 - 48.0938*x(i)*y(j) + 9.1571*x(i) + 2347.4087*y(j)^4 - 2066.9732*y(j)^3 + 498.3158*y(j)^2 - 34.3958*y(j) + 1.0559;


    end
end

[min_f, idx_f] = min(f(:));
fprintf('min_f3 = %g\n', min_f);


p = polyfitn([xvec',yvec'],zvec,5);
polyn2sympoly(p);
fprintf('RMSE - polynomial of order 5: %g\n', p.RMSE);

f = zeros(length(psi_range), length(phi_range));
for i = 1:length(x)
    for j = 1:length(y)
        f(i, j) = 4679.4338*x(i)^5 - 10034.2241*x(i)^4*y(j) - 4118.1666*x(i)^4 + 14674.7452*x(i)^3*y(j)^2 + 4807.141*x(i)^3*y(j) + 1252.0811*x(i)^3 - 11310.7334*x(i)^2*y(j)^3 - 6335.2789*x(i)^2*y(j)^2 - 50.8168*x(i)^2*y(j) - 168.1388*x(i)^2 - 41363.4115*x(i)*y(j)^4 + 24791.9896*x(i)*y(j)^3 - 2995.8425*x(i)*y(j)^2 + 53.4891*x(i)*y(j) + 11.5305*x(i) + 18676.7259*y(j)^5 + 1281.728*y(j)^4 - 4018.5308*y(j)^3 + 843.0939*y(j)^2 - 51.0421*y(j) + 1.2011;



    end
end

[min_f, idx_f] = min(f(:));
fprintf('min_f3 = %g\n', min_f);


p = polyfitn([xvec',yvec'],zvec,6);
polyn2sympoly(p);
fprintf('RMSE - polynomial of order 6: %g\n', p.RMSE);



f = zeros(length(psi_range), length(phi_range));
for i = 1:length(x)
    for j = 1:length(y)
        f(i, j) = 30443.6502*x(i)^6 + 31729.6485*x(i)^5*y(j) - 35025.9114*x(i)^5 + 5514.0011*x(i)^4*y(j)^2 - 42866.6728*x(i)^4*y(j) + 15664.849*x(i)^4 - 135132.2419*x(i)^3*y(j)^3 + 50803.2169*x(i)^3*y(j)^2 + 13716.5038*x(i)^3*y(j) - 3366.9682*x(i)^3 + 388032.6604*x(i)^2*y(j)^4 - 85444.4525*x(i)^2*y(j)^3 - 9638.6238*x(i)^2*y(j)^2 - 890.1771*x(i)^2*y(j) + 331.2147*x(i)^2 + 161927.7883*x(i)*y(j)^5 - 277540.3698*x(i)*y(j)^4 + 88327.9327*x(i)*y(j)^3 - 8247.7893*x(i)*y(j)^2 + 238.4126*x(i)*y(j) - 10.5159*x(i) - 752801.3339*y(j)^6 + 437971.9686*y(j)^5 - 74739.5959*y(j)^4 + 302.7975*y(j)^3 + 921.0363*y(j)^2 - 61.1967*y(j) + 1.4956;




    end
end

[min_f, idx_f] = min(f(:));
fprintf('min_f3 = %g\n', min_f);