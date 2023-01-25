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
pFile = strcat(diagpath, 'preassim_alp0.0_for1.0_T.nc');

% read state:
ncid = netcdf.open(pFile,'NC_NOWRITE');

[~, Ne] = netcdf.inqDim(ncid, 1);
[~, Nt] = netcdf.inqDim(ncid, 2);
[~, Nx] = netcdf.inqDim(ncid, 5);

xt   = squeeze(ncread(tFile, 'state'));
time = ncread(tFile, 'time'); 

forc = 1:15;     Nf = length(forc);
alp  = 0:.1:0.8; Na = length(alp);

for jj = 1:Na
    p1 = sprintf('%.1f', alp(jj));

    for ii = 1:Nf
        p2 = sprintf('%.1f', forc(ii)); 

        pFile = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_T.nc');
        aFile = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_T.nc');

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

        prior(jj, ii).mINF = mean(xim, 2);
        prior_mINF(jj, ii) = sum(prior(jj, ii).mINF) / Nx;

        analy(jj, ii).rmse = sqrt(mean( (xam - xt).^2 ));
        analy(jj, ii).sprd = sqrt(mean( (xas - 0 ).^2 ));
        analy(jj, ii).bias = mean( (xam - xt) );
        analy_mRMS(jj, ii) = sum(analy(jj, ii).rmse) / Nt;
        analy_mSTD(jj, ii) = sum(analy(jj, ii).sprd) / Nt;
        analy_mBAS(jj, ii) = sum(analy(jj, ii).bias) / Nt;

    end
end

for jj = 1:Na
    p1 = sprintf('%.1f', alp(jj));

    for ii = 1:Nf
        p2 = sprintf('%.1f', forc(ii)); 

        pFile = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_inf_T.nc');
        aFile = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_inf_T.nc');

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

        iprior(jj, ii).mINF = mean(xim, 2);
        iprior_mINF(jj, ii) = sum(iprior(jj, ii).mINF) / Nx;

        ianaly(jj, ii).rmse = sqrt(mean( (xam - xt).^2 ));
        ianaly(jj, ii).sprd = sqrt(mean( (xas - 0 ).^2 ));
        ianaly(jj, ii).bias = mean( (xam - xt) );
        ianaly_mRMS(jj, ii) = sum(ianaly(jj, ii).rmse) / Nt;
        ianaly_mSTD(jj, ii) = sum(ianaly(jj, ii).sprd) / Nt;
        ianaly_mBAS(jj, ii) = sum(ianaly(jj, ii).bias) / Nt;

    end
end

%%% for_every_obs: false

for jj = 1:Na
    p1 = sprintf('%.1f', alp(jj));

    for ii = 1:Nf
        p2 = sprintf('%.1f', forc(ii)); 

        pFile = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_F.nc');
        aFile = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_F.nc');

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

        fprior(jj, ii).rmse = sqrt(mean( (xfm - xt).^2 ));
        fprior(jj, ii).sprd = sqrt(mean( (xfs - 0 ).^2 ));
        fprior(jj, ii).bias = mean( (xfm - xt) );
        fprior_mRMS(jj, ii) = sum(fprior(jj, ii).rmse) / Nt;
        fprior_mSTD(jj, ii) = sum(fprior(jj, ii).sprd) / Nt;
        fprior_mBAS(jj, ii) = sum(fprior(jj, ii).bias) / Nt;

        fprior(jj, ii).mINF = mean(xim, 2);
        fprior_mINF(jj, ii) = sum(fprior(jj, ii).mINF) / Nx;

        fanaly(jj, ii).rmse = sqrt(mean( (xam - xt).^2 ));
        fanaly(jj, ii).sprd = sqrt(mean( (xas - 0 ).^2 ));
        fanaly(jj, ii).bias = mean( (xam - xt) );
        fanaly_mRMS(jj, ii) = sum(fanaly(jj, ii).rmse) / Nt;
        fanaly_mSTD(jj, ii) = sum(fanaly(jj, ii).sprd) / Nt;
        fanaly_mBAS(jj, ii) = sum(fanaly(jj, ii).bias) / Nt;

    end
end

for jj = 1:Na
    p1 = sprintf('%.1f', alp(jj));

    for ii = 1:Nf
        p2 = sprintf('%.1f', forc(ii)); 

        pFile = strcat(diagpath, 'preassim_alp', p1, '_for', p2, '_inf_F.nc');
        aFile = strcat(diagpath, 'analysis_alp', p1, '_for', p2, '_inf_F.nc');

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

        fiprior(jj, ii).rmse = sqrt(mean( (xfm - xt).^2 ));
        fiprior(jj, ii).sprd = sqrt(mean( (xfs - 0 ).^2 ));
        fiprior(jj, ii).bias = mean( (xfm - xt) );
        fiprior_mRMS(jj, ii) = sum(fiprior(jj, ii).rmse) / Nt;
        fiprior_mSTD(jj, ii) = sum(fiprior(jj, ii).sprd) / Nt;
        fiprior_mBAS(jj, ii) = sum(fiprior(jj, ii).bias) / Nt;

        fiprior(jj, ii).mINF = mean(xim, 2);
        fiprior_mINF(jj, ii) = sum(fiprior(jj, ii).mINF) / Nx;

        fianaly(jj, ii).rmse = sqrt(mean( (xam - xt).^2 ));
        fianaly(jj, ii).sprd = sqrt(mean( (xas - 0 ).^2 ));
        fianaly(jj, ii).bias = mean( (xam - xt) );
        fianaly_mRMS(jj, ii) = sum(fianaly(jj, ii).rmse) / Nt;
        fianaly_mSTD(jj, ii) = sum(fianaly(jj, ii).sprd) / Nt;
        fianaly_mBAS(jj, ii) = sum(fianaly(jj, ii).bias) / Nt;

    end
end


%% 

figure('pos', [220, 500, 850, 560])

ai = 6; %50%

subplot(2,2,[1, 3])
plot(forc, prior_mRMS(1, :), '-o', 'Color', bL, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto'); hold on; grid on 
plot(forc, prior_mRMS(ai, :), '-o', 'Color', rD, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto');

plot(forc, prior_mSTD(1, :), ':o', 'Color', bL, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto')
plot(forc, prior_mSTD(ai, :), ':o', 'Color', rD, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto')

yy = get(gca, 'YLim');
plot([8, 8], yy, '--', 'Color', gY, 'LineWidth', 1.5)

set(gca, 'FontSize', 16, 'XLim', [forc(1), forc(Nf)], 'XTick', 1:Nf, 'XTickLabel', forc, 'XGrid', 'off')
ylabel('Prior Diagnostics', 'FontSize', 18)
xlabel('Perturbed Forcing', 'FontSize', 18)

legend('EnKF RMSE', 'RD-EnKF RMSE', 'EnKF Spread', 'RD-EnKF Spread', 'Location', 'NorthWest', 'box', 'off')

title('No Inflation Case', 'FontSize', 22)

subplot(2,2,2)
plot(forc, iprior_mRMS(1, :), '-o', 'Color', bL, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto'); hold on; grid on 
plot(forc, iprior_mRMS(ai, :), '-o', 'Color', rD, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto');

plot(forc, iprior_mSTD(1, :), ':o', 'Color', bL, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto')
plot(forc, iprior_mSTD(ai, :), ':o', 'Color', rD, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto')

set(gca, 'FontSize', 16, 'XLim', [forc(1), forc(Nf)], 'XTick', 1:Nf, 'XTickLabel', forc, 'YLim', [0, 2], 'XGrid', 'off')
ylabel('Prior Diagnostics', 'FontSize', 18)

yy = get(gca, 'YLim');
plot([8, 8], yy, '--', 'Color', gY, 'LineWidth', 1.5)

title('Adaptive Inflation Case', 'FontSize', 22)

subplot(2,2,4)
plot(forc, iprior_mINF(1, :), '-o', 'Color', bL, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto'); hold on
plot(forc, iprior_mINF(ai, :), '-o', 'Color', rD, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'auto'); grid on

yy = get(gca, 'YLim');
plot([8, 8], yy, '--', 'Color', gY, 'LineWidth', 1.5)

set(gca, 'FontSize', 16, 'XLim', [forc(1), forc(Nf)], 'XTick', 1:Nf, 'XTickLabel', forc, 'XGrid', 'off')
ylabel('Average Prior Inflation', 'FontSize', 18)
xlabel('Perturbed Forcing', 'FontSize', 18)
