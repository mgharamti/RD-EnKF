clear 
clc

lB = [ 153, 255, 255 ]/255;
lP = [ 204, 153, 255 ]/255;
bL = [  30, 144, 255 ]/255;
rD = [ 255,  51,  51 ]/255;
gR = [   0, 153,   0 ]/255;
pR = [ 153,  51, 255 ]/255;
oR = [ 255, 153,  51 ]/255;
bK = [   0,   0,   0 ]/255;
gY = [ 150, 150, 150 ]/255;

figure('pos', [100, 200, 750, 360])

% p: prior p(x)
% l: likelihood p(y|x)
% q: posterior p(x|y)

y = 0;
m = 2;

sx = 1;
so = 0.5;

p1 = 1;
p2 = 1;

x = -3:.01:5;

a = 0.40;
c = 1/sqrt(2*pi);

p = c/sx * exp(-0.5 * (x-m).^2 / sx^2);
l = c/so * exp(-0.5 * (y-x).^2 / so^2);

k = a + (1-a) * l;

plot(x, p, 'Color', bL, 'LineWidth', 2); hold on
plot(x, l, 'Color', rD, 'LineWidth', 2); hold on
plot(x, k, 'Color', gR, 'LineWidth', 2); grid on

qe = p .* l;
qe = qe / trapz(x, qe);

qr = p .* k;
qr = qr / trapz(x, qr);

plot(x, qe, 'Color', gY, 'LineWidth', 3)
plot(x, qr, 'Color', bK, 'LineWidth', 3)

set(gca, 'FontSize', 16, 'YLim', [0, 1], 'YTick', 0:.1:1)
xlabel('Variable of interest', 'FontSize', 18)
ylabel('Distribution', 'FontSize', 18)

L = legend('Prior', 'Likelihood', 'Mixture Likelihood', 'EnKF Posterior', 'RD-EnKF Posterior', ...
    'Location', 'NorthWest');
title(L, 'Dormancy Rate: 0.1');

axes('Position', [.605, .563, .30, .27])

A = 0.02:0.02:0.98;
N = length(A);

Cols = cbrewer('qual', 'Accent', N, 'PCHIP'); box on
colormap(Cols); hold on
for i = 1:N
    pyx = A(i) + (1-A(i)) * l;
    pxy = p .* pyx;
    pxy = pxy / trapz(x, pxy);
    plot(x, pxy, 'Color', Cols(i, :))
end
set(gca, 'FontSize', 13, 'YLim', [0, 0.8], 'XLim', [-3, 5])
title('Possible RD-EnKF Posteriors', 'FontSize', 15)
h = colorbar;

set(h, 'YTick', A(1:8:end), 'Position', [.91, .565, .02, .265], 'FontSize', 12)
set(get(h, 'title'), 'String', '\alpha', 'FontSize', 14)

xx = [0.63 0.68];   % adjust length and location of arrow 
yy = [0.34 0.54];
annotation('textarrow', xx, yy, 'FontSize', 13, 'Linewidth', 1)

