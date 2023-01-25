clear 
clc

bK = [   0,   0,   0 ]/255;
bL = [  30, 144, 255 ]/255;
rD = [ 255,  51,  51 ]/255;
gR = [   0, 153,   0 ]/255;
pR = [ 153,  51, 255 ]/255;
oR = [ 255, 153,  51 ]/255;
gY = [.8, .8, .8];

volume = '/Volumes/DRIVE/Work/rd_enkf/l96/';
filter = 'EnKF_all/';
diags  = 'imperfect/';
diagpath = strcat(volume, filter, diags);

% diag files:
tFile = strcat(diagpath, '../true_state.nc');
pFile = strcat(diagpath, 'preassim_alp0.0_for8.0_T.nc');

% read state:
ncid = netcdf.open(pFile,'NC_NOWRITE');

[~, Ne] = netcdf.inqDim(ncid, 1);
[~, Nt] = netcdf.inqDim(ncid, 2);
[~, Nx] = netcdf.inqDim(ncid, 5);

xt   = squeeze(ncread(tFile, 'state'));
time = ncread(tFile, 'time'); 

forc = 1:15;     Nf = length(forc);
alp = 0:0.1:0.8; Na = length(alp);

F = 8;
for jj = 1:Na

    p1 = sprintf('%.1f', alp(jj));
    p2 = sprintf('%.1f', F);  

    pFile = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_T.nc');
    aFile = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_T.nc');
    pFili = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_inf_T.nc');
    aFili = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_inf_T.nc');

    xf(:, :, :, jj) = ncread(pFile, 'state');
    xa(:, :, :, jj) = ncread(aFile, 'state');
    yf(:, :, :, jj) = ncread(pFili, 'state');
    ya(:, :, :, jj) = ncread(aFili, 'state');

end

F = 1;
for jj = 1:Na

    p1 = sprintf('%.1f', alp(jj));
    p2 = sprintf('%.1f', F);  

    pFile = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_T.nc');
    aFile = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_T.nc');
    pFili = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_inf_T.nc');
    aFili = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_inf_T.nc');

    xf_1(:, :, :, jj) = ncread(pFile, 'state');
    xa_1(:, :, :, jj) = ncread(aFile, 'state');
    yf_1(:, :, :, jj) = ncread(pFili, 'state');
    ya_1(:, :, :, jj) = ncread(aFili, 'state');

end

F = 15;
for jj = 1:Na

    p1 = sprintf('%.1f', alp(jj));
    p2 = sprintf('%.1f', F);  

    pFile = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_T.nc');
    aFile = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_T.nc');
    pFili = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_inf_T.nc');
    aFili = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_inf_T.nc');

    xf_15(:, :, :, jj) = ncread(pFile, 'state');
    xa_15(:, :, :, jj) = ncread(aFile, 'state');
    yf_15(:, :, :, jj) = ncread(pFili, 'state');
    ya_15(:, :, :, jj) = ncread(aFili, 'state');

end

%%
x  = 20;
w  = 1;
h  = 6; 
yl = 0.2;

figure('pos', [2000, 500, 900, 900])

E_8 = rh(squeeze(xf(x, :, :, 1))', xt(x, :));
R_8 = rh(squeeze(xf(x, :, :, h))', xt(x, :));

E_1 = rh(squeeze(xf_1(x, :, :, 1))', xt(x, :));
R_1 = rh(squeeze(xf_1(x, :, :, h))', xt(x, :));

E_15 = rh(squeeze(xf_15(x, :, :, 1))', xt(x, :));
R_15 = rh(squeeze(xf_15(x, :, :, h))', xt(x, :));

subplot(321)
bar(1:21, E_1, w, 'FaceColor', bL, 'EdgeColor', 'k'); grid on

set(gca, 'FontSize', 14, 'XLim', [0.5, 21.5], 'XTick', 1:21, 'YLim', [0, yl], 'YTickLabelRotation', 90)
title('EnKF', 'FontSize', 20)
xlabel('Rank', 'FontSize', 18)
ylabel('Probability', 'FontSize', 18)
L = legend('F = 1', 'Location', 'North');
set(L, 'FontSize', 16)

subplot(322)
bar(1:21, R_1, w, 'FaceColor', bL, 'EdgeColor', 'k'); grid on

set(gca, 'FontSize', 14, 'XLim', [0.5, 21.5], 'XTick', 1:21, 'YLim', [0, yl], 'YTickLabelRotation', 90)
title(['RD-EnKF (' num2str(alp(h)*100) '%)'], 'FontSize', 20)
xlabel('Rank', 'FontSize', 18)
ylabel('Probability', 'FontSize', 18)
L = legend('F = 1', 'Location', 'North');
set(L, 'FontSize', 16)

subplot(323)
bar(1:21, E_8, w, 'FaceColor', bL, 'EdgeColor', 'k'); grid on

set(gca, 'FontSize', 14, 'XLim', [0.5, 21.5], 'XTick', 1:21, 'YLim', [0, yl], 'YTickLabelRotation', 90)
xlabel('Rank', 'FontSize', 18)
ylabel('Probability', 'FontSize', 18)
L = legend('F = 8', 'Location', 'North');
set(L, 'FontSize', 16)

subplot(324)
bar(1:21, R_8, w, 'FaceColor', bL, 'EdgeColor', 'k'); grid on

set(gca, 'FontSize', 14, 'XLim', [0.5, 21.5], 'XTick', 1:21, 'YLim', [0, yl], 'YTickLabelRotation', 90)
xlabel('Rank', 'FontSize', 18)
ylabel('Probability', 'FontSize', 18)
L = legend('F = 8', 'Location', 'North');
set(L, 'FontSize', 16)

subplot(325)
bar(1:21, E_15, w, 'FaceColor', bL, 'EdgeColor', 'k'); grid on

set(gca, 'FontSize', 14, 'XLim', [0.5, 21.5], 'XTick', 1:21, 'YLim', [0, yl], 'YTickLabelRotation', 90)
xlabel('Rank', 'FontSize', 18)
ylabel('Probability', 'FontSize', 18)
L = legend('F = 15', 'Location', 'North');
set(L, 'FontSize', 16)

subplot(326)
bar(1:21, R_15, w, 'FaceColor', bL, 'EdgeColor', 'k'); grid on

set(gca, 'FontSize', 14, 'XLim', [0.5, 21.5], 'XTick', 1:21, 'YLim', [0, yl], 'YTickLabelRotation', 90)
xlabel('Rank', 'FontSize', 18)
ylabel('Probability', 'FontSize', 18)
L = legend('F = 15', 'Location', 'North');
set(L, 'FontSize', 16)
