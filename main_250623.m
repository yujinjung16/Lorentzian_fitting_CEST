% Script for Lorentzian fitting on CEST data acquired using amide phantom
% Contact : Yujin Jung <emily7962@kaist.ac.kr>

%% Prepare
clear; clc; close all;
dirHome = 'E:\Lorentzian_fitting_250623';
load(fullfile(dirHome, 'data'),'xdata','ydata')
pools = {'Water','MT','Amide','2.5 ppm','4.5 ppm','-4.5 ppm'};

%% Fitting
%             1 = water;        2 = MT;           3 = amide;        4 = ? (2.5ppm);   5 = ? (4.5ppm);   6 = NOE? (-4.5ppm);
%       z1    A1    F1    w1    A2    F2    w2    A4    F4    w4    A5    F5    w5    A6    F6    w6    A7    F7    w7
lb = [  0     0     0     -0.5  0     10    -3    0     0     3.0   0     0     2.0   0     0     4.0   0     0     -5.0];
x0 = [  0.95  1     0.2   0     0.05  25    -2    0.07  1.0   3.3   0.09  0.9   2.5   0.015 0.5   4.5   0.02  0.5   -4.5];
ub = [  1     1     1.5   0.5   0.2   100   -1    0.2   1.5   4.0   0.2   1.5   3.0   0.05  1.0   5.0   0.05  1.0   -4.0];

fun = @(p,xdata)p(1)  -  p(2) .* (p(3).^2./4) ./ ((xdata-p(4)).^2 + (p(3).^2./4))... % pool 1
    -  p(5) .* (p(6).^2./4) ./ ((xdata-p(7)).^2 + (p(6).^2./4))... % pool 14
    -  p(8) .* (p(9).^2./4) ./ ((xdata-p(10)).^2 + (p(9).^2./4))... % pool 3
    -  p(11) .* (p(12).^2./4) ./ ((xdata-p(13)).^2 + (p(12).^2./4))... % pool 4
    -  p(14) .* (p(15).^2./4) ./ ((xdata-p(16)).^2 + (p(15).^2./4))... % pool 5
    -  p(17) .* (p(18).^2./4) ./ ((xdata-p(19)).^2 + (p(18).^2./4)); % pool 6
fun_1pool = @(p,xdata)1  -  p(2) .* (p(3).^2./4) ./ ((xdata-p(4)).^2 + (p(3).^2./4));
coeffs = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);

%% Display z-spectrum w/ measurements and fitting results
figure(1);
xfit = linspace(xdata(2),xdata(end-1));
yfit = fun(coeffs, xfit);
hold on;
plot(xdata(2:50),ydata(2:50),'o','MarkerSize', 5,'color',[0 0.4470 0.7410],'Displayname','Meas')
plot(xfit,yfit,'LineWidth', 2, 'color',[0 0.4470 0.7410],'Displayname','Fit')

set(gca,'Xdir','reverse')
xline(3.5,'--r','HandleVisibility','off')
ylim([0, 1.05])
grid on;
grid minor;

ax = gca;
ax.FontSize = 11;
xlabel('Frequency offset [ppm]','FontSize',13);
ylabel('S/S_0')
legend('Location','SouthEast','FontSize',11)
title(sprintf('Z-spectrum'),'FontSize',14);

%% Display each pool
f = figure(2);s
xfit = linspace(xdata(2),xdata(end-1));
yfit = fun(coeffs, xfit);
subplot(2,4,1); hold on;
plot(xdata(2:50),ydata(2:50),'.','MarkerSize', 3,'color',[0 0.4470 0.7410],'HandleVisibility','off')
plot(xfit,yfit,'LineWidth', 2, 'color',[0 0.4470 0.7410])

for n1 = 1 : length(pools)
    coeffs_1pool = [coeffs(1), coeffs((n1-1)*3+2:(n1-1)*3+4)];
    yfit = fun_1pool(coeffs_1pool, xfit);
    subplot(2,4,n1+1); hold on;
    plot(xfit,yfit,'LineWidth', 2, 'color',[0 0.4470 0.7410])
end

for n2 = 1 : length(pools)+1
    subplot(2,4,n2);    
    set(gca,'Xdir','reverse')
    grid on;
    grid minor;
    xlabel('Frequency offset [ppm]');
    
    if n2 == 1
        xline(3.5,'--r','HandleVisibility','off')
        ylim([0, 1.05])
        ylabel('S/S_0');
        title(sprintf('Z-spectrum'));
    else
        pool = pools{n2-1};
        ylabel('Fitted');
        title(pool);
    end    
end

f.Position = [200 500 1100 500];