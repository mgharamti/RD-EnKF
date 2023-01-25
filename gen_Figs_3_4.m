clear 

rng(2)

gY = [ 150, 150, 150 ]/255;
lB = [ 153, 255, 255 ]/255;
lP = [ 204, 153, 255 ]/255;
bK = [   0,   0,   0 ]/255;
bL = [  30, 144, 255 ]/255;
rD = [ 255,  51,  51 ]/255;
gR = [   0, 153,   0 ]/255;
pR = [ 153,  51, 255 ]/255;
oR = [ 255, 153,  51 ]/255;

addpath ../

nonlinear = false;
if nonlinear 
    params = 2; %#ok
else
    params = 1; %#ok
end

N = 10; 
T = 20; %100
S = 1:N;

E  = randn(1, N);

A = [0.0, -0.3];
B = [2, 2];

a = A(params);
b = B(params);
r = 0.1;

Y = randn(1, T);

%% 
alpha_cases = [0, 0.1];

figure('pos', [2500, 430, 900, 600])

tiledlayout(2,1)

for i = 1:2
    alpha = alpha_cases(i); 

    Nd = round(alpha * N);       
    Na = N - Nd;  

    Ea = E;
    
    Xf = zeros(N, T);
    Xa = zeros(N, T+1);

    Xa(:, 1) = Ea';
    for t = 1:T
    
        % Now advance 
        Ef = move_forward(Ea, a, b);
        Ea = Ef;
    
        % observation 
        y = Y(t);  
        
        % randomization
        active_mem = sort(randperm(N, Na));
        Er = Ea(active_mem);
    
        % update
        dx = obs_increment_eakf(Er, y, r);
        Er = Er + dx;
    
        % reunion
        Ea(active_mem) = Er; 

        % Save for diags
        Xf(:, t)   = Ef';
        Xa(:, t+1) = Ea';
        Da(t, :)   = S( ismember(1:N, active_mem));
        Dm(t, :)   = S(~ismember(1:N, active_mem));
    end

    REf(:, :, i) = Xf;
    REa(:, :, i) = Xa;

    nexttile

    hold on 
    grid on

    step = 0.1;
    ax1  = gca;
    
    truth = plot([0, T], [0, 0], '--', 'Color', 'k', 'LineWidth', 2);
    
    plot(0, 0, '*', 'Color', rD, 'MarkerSize', 12)
    plot(step, Xa(:, 1), '*', 'Color', bL, 'MarkerSize', 6)
    
    set(gca, 'FontSize', 18, 'XLim', [0, T], 'XTick', 1:T-1, 'YLim', [-10, 8], 'box', 'on')
    
    if i > 1, xlabel('Time Steps', 'FontSize', 20); end
    ylabel('State', 'FontSize', 20)
    
    for t = 1:T
        lx = (t - 1) + step;
        rx = t - step;
        ox = t;
        nx = t + step;
    
        bx(1, 1:N) = lx;
        bx(2, 1:N) = rx;
    
        by(1, 1:N) = Xa(:, t);
        by(2, 1:N) = Xf(:, t);
    
        if Nd > 0
            pr=plot(bx(:, Da(t, :)), by(:, Da(t, :)), '-', 'Color', gR);
            pd=plot(bx(:, Dm(t, :)), by(:, Dm(t, :)), '-', 'Color', oR);
            plot(bx(2, Da(t, :)), by(2, Da(t, :)), '*', 'Color', gR)
            plot(bx(2, Dm(t, :)), by(2, Dm(t, :)), '*', 'Color', oR)
        else
    
            plot(bx, by, 'Color', gR)
            pr = plot(bx(2, :), by(2, :), '*', 'Color', gR);
        end
    
        ob = plot(t, Y(t), '*', 'Color', rD, 'MarkerSize', 12);
    
        po = plot(nx, Xa(:, t+1), '*', 'Color', bL, 'MarkerSize', 6);
    end
    
    text(2.05, -6, sprintf('Prior Error: %.3f, Posterior Error: %.3f', mean(abs(mean(Xf, 1) - Y)), mean(abs(mean(Xa, 1) - [0, Y]))), ...
            'FontSize', 18)
    
    if alpha == 0, msg = 'EnKF'; end
    if alpha ~= 0, msg = sprintf('RD-EnKF'); end
    
    title(msg, 'FontSize', 24, 'FontWeight', 'Bold')
    
    yyaxis right
    ss = std(Xa, 0, 1); 
    ps = plot(0:t, ss, '-o', 'LineWidth', 1.5, 'Color', gY, 'MarkerSize', 8, 'MarkerFaceColor', 'auto');
    set(gca, 'YColor', gY, 'YScale', 'log', 'YLim', [1.e-3, 5], 'YTick', [1.e-3, 1.e-2, 1e-1, 1])
    ylabel('Spread', 'FontSize', 20)

    text(3.75, 0.0023, sprintf('Ensemble Spread: %.3f', mean(ss)), 'FontSize', 18)

    clear Da Dm
end

%% Table 1: 
tt = [1, 2, 3, 4, 5, 10, 20]; %, 50, 100, 1000];

Xef = REf(:, :, 1);
Xrf = REf(:, :, 2);

Xea = REa(:, :, 1);
Xra = REa(:, :, 2);

Z = [0, Y];
for i = 1:length(tt)

    eef = mean(abs(mean(Xef(:, 1:tt(i)), 1) - Y(1:tt(i))));
    rrf = mean(abs(mean(Xrf(:, 1:tt(i)), 1) - Y(1:tt(i))));

    eea = mean(abs(mean(Xea(:, 1:tt(i)), 1) - Z(1:tt(i))));
    rra = mean(abs(mean(Xra(:, 1:tt(i)), 1) - Z(1:tt(i))));

    disp(['After ' num2str(tt(i)) ' cycle, f_imp is ' num2str(100 - rrf/eef * 100) '%, ' ...
          'a_imp is ' num2str(100 - rra/eea * 100) '%'])

end
disp(['After ' num2str(1000) ' cycle, f_imp is ' num2str(100 - 1.929/2.961 * 100) '%, ' ...
      'a_imp is ' num2str(100 -0.518/1.315 * 100) '%'])
