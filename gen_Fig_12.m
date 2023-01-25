clear 
clc

gY = [ 150, 150, 150 ]/255;

obs  = [2231342, 2231600, 2232155, 2232000, 2232200, 2232400, 2232500, 2233475, ...
        2233484, 2233500, 2234000, 2234324, 2234344, 2234384, 2234400, 2234990, ...
        2235000, 2235200, 2235500, 2236500, 2236605, 2237293, 2238500, 2239000, ...
        2239501, ...
        ...
        2240000, 2240500, 2243000, 2243960, 2244333, 2245500, 2247510, 2248000, ...
        2249007, 2249500, ...
        ...
        2251000, 2251500, 2256500, ...
        ...
        2262900, 2263800, 2263869, 2264003, 2264030, 2264051, 2264100, 2264495, ...
        2266200, 2266291, 2266300, 2266480, 2266496, ...
        ...
        2270000, 2270500, 2271500, 2272650, 2272676, 2274005, 2274010, 2274490, ...
        2275197, ...
        ...
        2291580, 2292010, 2293987, 2294161, 2294217, 2294330, 2294405, 2294650, ...
        2294898, 2295194, 2295420, 2295520, 2295580, 2295637, 2296260, 2296500, ...
        2296750, 2297100, 2297155, 2297310, 2297600, 2298110, 2298123, 2298202, ...
        2298488, 2298492, 2298495, 2298527, 2298530, 2298554, 2298608, 2298760, ...
        2298830, 2298880, 2299410, 2299450, 2299472, 2299950, ...
        ...
        2300017, 2300018, 2300033, 2300042, 2300075, 2300100, 2300210, 2300300, ...
        2300500, 2300700, 2300882, 2301000, 2301150, 2301300, 2301500, 2301738, ...
        2301740, 2301745, 2301750, 2301793, 2301900, 2301990, 2302500, 2303000, ...
        2303330, 2303350, 2304500, 2306500, 2306647, 2306774, 2306950, 2307000, ...
        2307323, 2307359, 2307445, 2307668, 2307674, 2307780, 2309415, 2309421, ...
        2309425, 2309740, ...
        ...
        2310000, 2310525, 2310675, 2310947, 2311000, 2311500, 2312000, 2312180, ...
        2312200, 2312300, 2312500, 2312558, 2312598, 2312600, 2312640, 2312645, ...
        2312667, 2312675, 2312700, 2312720, 2312762, 2312764, 2313000, 2313098, ...
        2313100, ...
        ...
        2320700, 2321000, 2321500, 2321898, 2321958, 2322700, 2322800, ...
        ...
        ];

dir_ol   = {'Exp/N80_ol'};
dir_exps = {'Exp/N20_infprpo', ...
            'Exp_rd/rd_N20_infprpo_a0.1', 'Exp_rd/rd_N20_infprpo_a0.2', ...
            'Exp_rd/rd_N20_infprpo_a0.3', 'Exp_rd/rd_N20_infprpo_a0.4', ...
            'Exp_rd/rd_N20_infprpo_a0.5'};

[~, ~, ~, ~, obs_all, openloop_all, forecast_all, analysis_all] = HydroDARTdiags(dir_exps, obs, dir_ol, 0, 0, 0);

% separate exp name from path
exp_name = string(missing);
for e = 1:length(dir_exps)
    sp_names    = strsplit(dir_exps{e}, '/');
    exp_name(e) = sp_names{end};
end

exp_names_labels = ['OL', exp_name];

%% Compute RMSE-based on percentiles

any_exp = 1;
obs_val = 6;
metrics = 2;
bias_dg = 3;

pr = 10; %10th percentile
Ny = size(forecast_all, 2); 
Ne = length(dir_exps);

OBS_p = zeros(Ny, 1);
for e = 1:Ne+1
    for k = 1:Ny  
        gage_dat = obs_all{obs_val, k, any_exp};
        OBS_p(k) = prctile(gage_dat, pr);
    
        low_flow = gage_dat <= OBS_p(k); 
        hig_flow = gage_dat > OBS_p(k);
    
        if e == 1
            rmse.f = openloop_all{metrics, k};
            rmse.a = rmse.f;
        else
            rmse.f = forecast_all{metrics, k, e-1};
            rmse.a = analysis_all{metrics, k, e-1};
            bias.f = forecast_all{bias_dg, k, e-1};
            bias.a = analysis_all{bias_dg, k, e-1};
        end
    
        p_RMSE_LOW.f(k, e) = mean(rmse.f(low_flow));
        p_RMSE_LOW.a(k, e) = mean(rmse.a(low_flow));

        p_RMSE_HIG.f(k, e) = mean(rmse.f(hig_flow));
        p_RMSE_HIG.a(k, e) = mean(rmse.a(hig_flow));

        p_RMSE_ALL.f(k, e) = mean(rmse.f);
        p_RMSE_ALL.a(k, e) = mean(rmse.a);

        if e > 1
            p_BIAS_ALL.f(k, e-1) = mean(abs(bias.f));
            p_BIAS_ALL.a(k, e-1) = mean(abs(bias.a));
        end
    end
end

ens_val = 8;
ens_avg = 6;
ole_val = 4;
spr_val = 4;

% The function crps can be downloaded from:
% https://www.mathworks.com/matlabcentral/fileexchange/47807-continuous-rank-probability-score

CRPSS = zeros(Ne, Ny); 
Esprf = zeros(Ne, Ny); 
Espra = Esprf;
for k = 1:Ny
    obs = obs_all{obs_val, k, 1};

    for e = 1:Ne
        pr_ens = forecast_all{ens_val, k, e};
        ol_ens = openloop_all{ole_val, k};

        exp_CRPS = crps(pr_ens', obs);
        ref_CRPS = crps(ol_ens', obs);

        CRPSS(e, k) = 1 - exp_CRPS / ref_CRPS;

        Esprf(e, k) = nanmean(forecast_all{spr_val, k, e});
        Espra(e, k) = nanmean(analysis_all{spr_val, k, e});
    end
end

%% CRPSS 
figure('pos', [285, 420, 550, 520]) 

Ns = size(CRPSS, 1);
CRPSS(CRPSS<0) = 0;

COLS = turbo(Ne+1);

boxplot(CRPSS', 'Notch', 'on', 'Whisker', 1.5, 'Widths', 0.5, 'OutlierSize', 12, ...
    'Colors', COLS(2:end, :), 'Symbol', '*')

set(gca, 'FontSize', 16, 'XLim', [.5, Ns+0.5], 'XTick', 1:Ns, 'XTickLabel', ...
    {'EnKF', 'RD-EnKF, 10%', 'RD-EnKF, 20%', 'RD-EnKF, 30%', ...
     'RD-EnKF, 40%', 'RD-EnKF, 50%'}, 'XGrid', 'on', 'YLim', [-0.15, 1], ...
     'YTick', 0:.2:1, 'YTickLabelRotation', 90)

ylabel('Continuous Ranked Probability Skill Score', 'FontSize', 18)

mCRPSS = mean(CRPSS, 2);

yy = -0.1; offset = 0.17;
for k = 1:Ne
    text(k-offset, yy, sprintf('%.2f', mCRPSS(k)), 'FontSize', 16, 'Color', COLS(k+1, :), 'FontWeight', 'bold');
end

yyaxis right

plot(1:Ne, mean(Esprf, 2), '-o', 'LineWidth', 3, 'Color', gY, ... 
     'MarkerSize', 12, 'MarkerFaceColor', 'auto'); hold on
plot(1:Ne, mean(Espra, 2), '--o', 'LineWidth', 3, 'Color', gY, ... 
     'MarkerSize', 12, 'MarkerFaceColor', 'auto')

set(gca, 'FontSize', 16, 'YColor', gY, 'YLim', [-1, 6], 'YTick', 0:1:5)
ylabel('Ensemble Spread', 'FontSize', 18)
