clear all;
clc;
load variables\exsearch100\paramsW.mat

fontsize = 16;
% plot error vs. psi and phi
fig = figure();clf;
subplot(1, 2, 1)
surf(psi_range, phi_range, error);
colormap(jet);
xlabel('\psi', 'fontsize',fontsize);
ylabel('\phi', 'fontsize',fontsize);
zlabel('NRMSE', 'fontsize',fontsize);
title('Error vs. optimization parameters for W', 'fontsize',16);
hold on;
set(gca,'fontsize',fontsize);
grid on;
subplot(1, 2, 2)
% contour plot
contourf(psi_range, phi_range, error);
colormap(jet);
cbr = colorbar('location','southoutside');
set(cbr,'ZTick',0:0.1:1.5)
xlabel('\psi', 'fontsize',fontsize);
ylabel('\phi', 'fontsize',fontsize);
title('Error contour plot for W', 'fontsize',16);
hold on;
set(gca,'fontsize',fontsize);
grid on;
saveas(fig, strcat('images/surfContW.png'), 'png');
saveas(fig, strcat('images/surfContW.fig'), 'fig');


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
        f(i, j) = 1.3477.*x(i).^3 + 36.7252.*x(i).^2.*y(j) - 2.7548.*x(i).^2 - 173.8875.*x(i).*y(j).^2 + 30.5263.*x(i).*y(j) - 0.030841.*x(i) + 519.6074.*y(j).^3 - 119.2905.*y(j).^2 - 1.3545.*y(j) + 1.4845;
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
title('Error approximation contour plot for W', 'fontsize',16);
hold on;
set(gca,'fontsize',fontsize);
grid on;
saveas(fig, strcat('images/Waprox3.png'), 'png');
saveas(fig, strcat('images/Waprox3.fig'), 'fig');


p = polyfitn([xvec',yvec'],zvec,4);
polyn2sympoly(p);
fprintf('RMSE - polynomial of order 4: %g\n', p.RMSE);
for i = 1:length(x)
    for j = 1:length(y)
        f(i, j) = 98.086*x(i)^4 + 111.6503*x(i)^3*y(j) - 88.2862*x(i)^3 + 60.6086*x(i)^2*y(j)^2 - 42.3867*x(i)^2*y(j) + 24.4551*x(i)^2 - 1551.3303*x(i)*y(j)^3 + 267.2682*x(i)*y(j)^2 + 8.9945*x(i)*y(j) - 2.4223*x(i) - 1765.6803*y(j)^4 + 1536.1456*y(j)^3 - 301.2762*y(j)^2 + 9.3638*y(j) + 1.4051;
    end
end
[min_f, idx_f] = min(f(:));
fprintf('min_f3 = %g\n', min_f);
[min_err, idx_error] = min(error(:));
[min_f, idx_f] = min(f(:));
err_aprox = error(m, n);
rel_err = (abs(min_err - err_aprox)/min_err)*100;
fprintf('Relative error (in percentage) of polyn. of order 4\n to the exhaustive search minimum: %g \n', rel_err );

figure();clf;
surf(x, y, f);
title('Polynome of order 4', 'fontsize',16);
xlabel('\psi', 'fontsize',16);
ylabel('\phi', 'fontsize',16);
zlabel('NRMSE', 'fontsize',16);


p = polyfitn([xvec',yvec'],zvec,5);
polyn2sympoly(p)
fprintf('RMSE - polynomial of order 5: %g\n', p.RMSE);
for i = 1:length(x)
    for j = 1:length(y)
        f(i, j) = -418.6279*x(i)^5 + 2345.3541*x(i)^4*y(j) + 282.1785*x(i)^4 + 2772.9737*x(i)^3*y(j)^2 - 2319.2277*x(i)^3*y(j) - 30.8391*x(i)^3 - 9853.9594*x(i)^2*y(j)^3 + 1353.0122*x(i)^2*y(j)^2 + 535.9372*x(i)^2*y(j) - 8.6812*x(i)^2 + 43883.8993*x(i)*y(j)^4 - 15163.3063*x(i)*y(j)^3 + 1599.0783*x(i)*y(j)^2 - 90.9155*x(i)*y(j) + 1.7263*x(i) - 73543.2482*y(j)^5 + 26229.1639*y(j)^4 - 1734.4303*y(j)^3 - 196.2443*y(j)^2 + 11.8456*y(j) + 1.2795;
    end
end
[min_f, idx_f] = min(f(:));
fprintf('min_f3 = %g\n', min_f);

p = polyfitn([xvec',yvec'],zvec,6);
polyn2sympoly(p);
fprintf('RMSE - polynomial of order 6: %g\n', p.RMSE);
for i = 1:length(x)
    for j = 1:length(y)
        f(i, j) = -14213.8134.*x(i).^6 - 13250.6662.*x(i).^5.*y(j) + 17963.0148.*x(i).^5 + 21856.0335.*x(i).^4.*y(j).^2 + 11224.8136.*x(i).^4.*y(j) - 8636.7178.*x(i).^4 - 3074.9576.*x(i).^3.*y(j).^3 - 13789.3657.*x(i).^3.*y(j).^2 - 3595.5023.*x(i).^3.*y(j) + 1966.9785.*x(i).^3 - 27972.0113.*x(i).^2.*y(j).^4 + 3179.8197.*x(i).^2.*y(j).^3 + 3846.9361.*x(i).^2.*y(j).^2 + 446.5152.*x(i).^2.*y(j) - 212.619.*x(i).^2 + 246233.3169.*x(i).*y(j).^5 - 68043.9546.*x(i).*y(j).^4 + 1901.2664.*x(i).*y(j).^3 + 241.5359.*x(i).*y(j).^2 - 38.0555.*x(i).*y(j) + 9.4186.*x(i) + 708299.2404.*y(j).^6 - 547769.4558.*y(j).^5 + 146514.4439.*y(j).^4 - 16024.4582.*y(j).^3 + 600.8752.*y(j).^2 - 6.1947.*y(j) + 1.3054;
    end
end

[min_f, idx_f] = min(f(:));
fprintf('min_f3 = %g\n', min_f);
[min_err, idx_error] = min(error(:));
[min_f, idx_f] = min(f(:));
[m, n] = ind2sub(size(f),idx_f);
err_aprox = error(m, n);
rel_err = (abs(min_err - err_aprox)/min_err)*100;
fprintf('Relative error (in percentage) of polyn. of order 6\n to the exhaustive search minimum: %g \n', rel_err );

figure();clf;
surf(x, y, f);
title('Polynome of order 6', 'fontsize',16);
xlabel('\psi', 'fontsize',16);
ylabel('\phi', 'fontsize',16);
zlabel('NRMSE', 'fontsize',16);


