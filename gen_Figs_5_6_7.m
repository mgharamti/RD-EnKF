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
diags  = 'perfect/';
diagpath = strcat(volume, filter, diags);

% diag files:
tFile = strcat(diagpath, '../true_state.nc');
pFile = strcat(diagpath, 'preassim_alp0.00_loc0.10.nc');

% read state:
ncid = netcdf.open(pFile,'NC_NOWRITE');

[~, Ne] = netcdf.inqDim(ncid, 1);
[~, Nt] = netcdf.inqDim(ncid, 2);
[~, Nx] = netcdf.inqDim(ncid, 5);

xt   = squeeze(ncread(tFile, 'state'));
time = ncread(tFile, 'time'); 

cut = .1:.02:0.3; Nl = length(cut);
alp = .0:.05:0.8; Na = length(alp);

% i- No inflation
for jj = 1:Na
    p1 = sprintf('%.2f', alp(jj));

    for ii = 1:Nl
        p2 = sprintf('%.2f', cut(ii)); 

        pFile = strcat(diagpath, 'preassim_alp', p1, '_loc', p2, '.nc');
        aFile = strcat(diagpath, 'analysis_alp', p1, '_loc', p2, '.nc');

        if ii == 1 && jj == 5
            xfe_rd4 = ncread(pFile, 'state');
            xae_rd4 = ncread(aFile, 'state');
        end

        xfm = ncread(pFile, 'state_mean');
        xfs = ncread(pFile, 'state_sd');
        xim = ncread(pFile, 'state_priorinf_mean');
        xam = ncread(aFile, 'state_mean');
        xas = ncread(aFile, 'state_sd');

        % If the run failed
        if size(xfm, 2) < size(xt, 2)
            xfm = NaN * ones(size(xt));
            xfs = xfm; 
            xim = xfm;
            xam = xfm;
            xas = xfm;
        end

        prior(jj, ii).rmse = sqrt(mean( (xfm - xt).^2 ));
        prior(jj, ii).sprd = sqrt(mean( (xfs - 0 ).^2 ));
        prior(jj, ii).bias = mean( (xfm - xt) );
        prior_mRMS(jj, ii) = sum(prior(jj, ii).rmse) / Nt;
        prior_mSTD(jj, ii) = sum(prior(jj, ii).sprd) / Nt;
        prior_mBAS(jj, ii) = sum(prior(jj, ii).bias) / Nt;

        prior(jj, ii).mINF = mean(xim, 1);

        analy(jj, ii).rmse = sqrt(mean( (xam - xt).^2 ));
        analy(jj, ii).sprd = sqrt(mean( (xas - 0 ).^2 ));
        analy(jj, ii).bias = mean( (xam - xt) );
        analy_mRMS(jj, ii) = sum(analy(jj, ii).rmse) / Nt;
        analy_mSTD(jj, ii) = sum(analy(jj, ii).sprd) / Nt;
        analy_mBAS(jj, ii) = sum(analy(jj, ii).bias) / Nt;

    end
end

% ii- Adding inflation
for jj = 1:Na
    p1 = sprintf('%.2f', alp(jj));

    for ii = 1:Nl
        p2 = sprintf('%.2f', cut(ii)); 

        pFile = strcat(diagpath, 'preassim_alp', p1, '_loc', p2, '_inf.nc');
        aFile = strcat(diagpath, 'analysis_alp', p1, '_loc', p2, '_inf.nc');

        xfm = ncread(pFile, 'state_mean');
        xfs = ncread(pFile, 'state_sd');
        xim = ncread(pFile, 'state_priorinf_mean');
        xam = ncread(aFile, 'state_mean');
        xas = ncread(aFile, 'state_sd');

        % If the run failed
        if size(xfm, 2) < size(xt, 2)
            xfm = NaN * ones(size(xt));
            xfs = xfm; 
            xim = xfm;
            xam = xfm;
            xas = xfm;
        end

        iprior(jj, ii).rmse = sqrt(mean( (xfm - xt).^2 ));
        iprior(jj, ii).sprd = sqrt(mean( (xfs - 0 ).^2 ));
        iprior(jj, ii).bias = mean( (xfm - xt) );
        iprior_mRMS(jj, ii) = sum(iprior(jj, ii).rmse) / Nt;
        iprior_mSTD(jj, ii) = sum(iprior(jj, ii).sprd) / Nt;
        iprior_mBAS(jj, ii) = sum(iprior(jj, ii).bias) / Nt;

        iprior(jj, ii).mINF = mean(xim, 1);

        ianaly(jj, ii).rmse = sqrt(mean( (xam - xt).^2 ));
        ianaly(jj, ii).sprd = sqrt(mean( (xas - 0 ).^2 ));
        ianaly(jj, ii).bias = mean( (xam - xt) );
        ianaly_mRMS(jj, ii) = sum(ianaly(jj, ii).rmse) / Nt;
        ianaly_mSTD(jj, ii) = sum(ianaly(jj, ii).sprd) / Nt;
        ianaly_mBAS(jj, ii) = sum(ianaly(jj, ii).bias) / Nt;

    end
end

%% compute diag:

figure('pos', [200, 400, 1220, 630])

ai = 5; %alpha=0.2
al = 1; %cutoff=0.1

subplot(221)

plot(time, prior(1, al).rmse, 'Color', bL, 'LineWidth', 2); hold on 
plot(time, prior(1, al).sprd, 'Color', rD, 'LineWidth', 2); grid on 
plot(time, prior(1, al).bias, 'Color', gR, 'LineWidth', 2); 

set(gca, 'FontSize', 18)
ylabel('Prior Diagnostics', 'FontSize', 20)

legend(sprintf('RMSE: %.2f'  , prior_mRMS(1, al)), ...
       sprintf('Spread: %.2f', prior_mSTD(1, al)), ...
       sprintf('Bias: %.3f'  , prior_mBAS(1, al)), 'Location', 'NorthWest')

title('EnKF (no inflation)', 'FontSize', 20, 'FontWeight', 'bold')

subplot(223)

plot(time, iprior(1, al).rmse, 'Color', bL, 'LineWidth', 2); hold on 
plot(time, iprior(1, al).sprd, 'Color', rD, 'LineWidth', 2); grid on 
plot(time, iprior(1, al).bias, 'Color', gR, 'LineWidth', 2); 

set(gca, 'FontSize', 18)
xlabel('Time (days)', 'FontSize', 20)
ylabel('Prior Diagnostics', 'FontSize', 20)

ax1 = gca;

title('EnKF (with inflation)', 'FontSize', 20, 'FontWeight', 'bold')

yyaxis right ; gY = [.6, .6, .6];
plot(time, iprior(1, al).mINF, 'Color', gY, 'LineWidth', 1)
set(gca, 'YColor', gY)

legend(ax1, sprintf('RMSE: %.2f', iprior_mRMS(1, al)), ...
       sprintf('Spread: %.2f'   , iprior_mSTD(1, al)), ...
       sprintf('Bias: %.3f'     , iprior_mBAS(1, al)), ...
       sprintf('Inflation: %.2f', mean(iprior(1, al).mINF)), 'Location', 'NorthWest')

subplot(222)

plot(time, prior(ai, al).rmse, 'Color', bL, 'LineWidth', 2); hold on 
plot(time, prior(ai, al).sprd, 'Color', rD, 'LineWidth', 2); grid on 
plot(time, prior(ai, al).bias, 'Color', gR, 'LineWidth', 2); 

set(gca, 'FontSize', 18)

legend(sprintf('RMSE: %.2f'  , prior_mRMS(ai, al)), ...
       sprintf('Spread: %.2f', prior_mSTD(ai, al)), ...
       sprintf('Bias: %.3f'  , prior_mBAS(ai, al)), 'Location', 'NorthWest')

title(['RD-EnKF (no inflation, ' num2str(alp(ai)*100) '%)'], 'FontSize', 20, 'FontWeight', 'bold')

subplot(224)

plot(time, iprior(ai, al).rmse, 'Color', bL, 'LineWidth', 2); hold on 
plot(time, iprior(ai, al).sprd, 'Color', rD, 'LineWidth', 2); grid on 
plot(time, iprior(ai, al).bias, 'Color', gR, 'LineWidth', 2); 

set(gca, 'FontSize', 18)
xlabel('Time (days)', 'FontSize', 20)
ax1 = gca;

title(['RD-EnKF (with inflation, ' num2str(alp(ai)*100) '%)'], 'FontSize', 20, 'FontWeight', 'bold')

yyaxis right
plot(time, iprior(ai, al).mINF, 'Color', gY, 'LineWidth', 1)
set(gca, 'YColor', gY)

legend(ax1, sprintf('RMSE: %.2f', iprior_mRMS(ai, al)), ...
       sprintf('Spread: %.2f'   , iprior_mSTD(ai, al)), ...
       sprintf('Bias: %.3f'     , iprior_mBAS(ai, al)), ...
       sprintf('Inflation: %.2f', mean(iprior(ai, al).mINF)), 'Location', 'NorthWest')


%% local sens
figure('pos', [2000, 500, 750, 400])

subplot('Position', [0.1, 0.2, 0.8, 0.7])

x_cut = 1:Na-1;
x_cut = x_cut([1, 3:2:end]);
Nk = length(x_cut);

alpp = alp(x_cut) * 100;
numm = round(Ne * alp(x_cut));

labelArray{1,1} = '#';
labelArray{2,1} = '%';

for k = 1:Nk
    labelArray{1,k+1} = num2str(numm(k));
    labelArray{2,k+1} = num2str(alpp(k));
end

nnl = Nl-3;

CM = turbo(nnl);

set(0, 'DefaultAxesColorOrder', CM); hold on; grid on

for ii = 1:nnl
    H = plot(alp(x_cut), prior_mRMS(x_cut, ii), '-*', 'LineWidth', 2, 'MarkerSize', 10);
end
ax = gca;

colormap(CM); 

plot([-1,alp], 2.7*ones(length(alp)+1), '--', 'Color', gY, 'LineWidth', 1.5)
 
set(ax, 'FontSize', 14, 'XLim', [-0.05, alp(x_cut(end))+0.05], 'XTick', alp(x_cut), ...
        'YGrid', 'off', 'box', 'on') 

ax.XTickLabel = '';

yl1 = get(ax, 'YLim');
if yl1(1) > 0
    yl2 = [yl1(1) - 0.25, yl1(2)];
else
    yl2 = yl1;
end

set(gca, 'YLim', yl2)

hc = colorbar;
set(hc, 'YTick', linspace(0, 1, nnl), 'YTickLabel', cut(1:nnl))
ylabel(hc, 'Localization Cutoff', 'FontSize', 16)

xx = alp(x_cut);

for i = 1:length(x_cut)+1
    if i == 1
        text(-alp(2), yl2(1)-0.05, sprintf('%s\n%s', labelArray{:, i}), ...
        'horizontalalignment', 'center', 'verticalalignment', 'top', 'FontSize', 14); 
    else
        text(xx(i-1), yl2(1)-0.05, sprintf('%s\n%s', labelArray{:, i}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'top', 'FontSize', 14);    
    end
end

ylabel('Average Prior RMSE', 'FontSize', 18)

text(0.22, -0.13 * ( yl2(2) - yl2(1) ), 'Dormant Members', 'FontSize', 18)


%% spurious correlations

for s = 1:Nx
    pairs = [s*ones(1, 40); 1:40];

    disp(['Var: ' num2str(s)])
    for p = 1:Np
        x1 = pairs(1, p);
        x2 = pairs(2, p);
        for k = 1:Nt
            tmpen = corrcoef(xfe(x1, :, k), xfe(x2, :, k));
            tmp80 = corrcoef(x80(x1, :, k), x80(x2, :, k));
            tmpr2 = corrcoef(xfe_rd4(x1, :, k), xfe_rd4(x2, :, k));
        
            En_corr(s, p, k) = tmpen(1, 2);
            En80_cr(s, p, k) = tmp80(1, 2);
            R2_corr(s, p, k) = tmpr2(1, 2);
        end
    end
end

m_En_corr = mean(En_corr, 3);
m_En80_cr = mean(En80_cr, 3);
m_R2_corr = mean(R2_corr, 3);

figure('pos', [130, 634, 1220, 350])

tiledlayout(1,3)

ax1 = nexttile;
imagesc(m_En80_cr); colormap(ax1, turbo); colorbar
set(gca, 'FontSize', 16, 'XTick', 5:5:35, 'YTick', 5:5:35)
xlabel('State Variables', 'FontSize', 18)
ylabel('State Variables', 'FontSize', 18)
title('80-mem Ensemble Correlations', 'FontSize', 20)


ax2 = nexttile;
imagesc(m_En_corr - m_En80_cr); colormap(ax2, redblue); caxis([-0.1, 0.1]); colorbar('Location', 'EastOutside');
set(gca, 'FontSize', 16, 'XTick', 5:5:35, 'YTick', [])
xlabel('State Variables', 'FontSize', 18)
title({'10-mem EnKF', 'Correlations Difference'}, 'FontSize', 18)

ax3 = nexttile;
imagesc(m_R2_corr - m_En80_cr); colormap(ax3, redblue); caxis([-0.1, 0.1]); hc = colorbar('Location', 'EastOutside');
set(gca, 'FontSize', 16, 'XTick', 5:5:35, 'YTick', [])
xlabel('State Variables', 'FontSize', 18)
title({'10-mem RD-EnKF (40%)', 'Correlations Difference'}, 'FontSize', 18)
