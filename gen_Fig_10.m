clear 
clc

ROUTE_FILE     = 'RouteLink.nc';

LAT            = double(ncread(ROUTE_FILE,'lat'))';
LON            = double(ncread(ROUTE_FILE,'lon'))';

gauges_manual  = gauges_2_indices_subset(ROUTE_FILE);
num_gauges_all = length(gauges_manual);

gauges_in      = gauges_manual(:, 2);

fromIndices    = double(ncread(ROUTE_FILE, 'fromIndices'));
fromIndsStart  = double(ncread(ROUTE_FILE, 'fromIndsStart'));
fromIndsEnd    = double(ncread(ROUTE_FILE, 'fromIndsEnd'));

[~, n_links]   = netcdf.inqDim(netcdf.open(ROUTE_FILE,'NC_NOWRITE'), 0);
no_up_links    = fromIndsStart == 0;
num_up_links   = fromIndsEnd - fromIndsStart + 1; num_up_links(no_up_links) = 0;

for k = 1:length(gauges_in)
    tmp = gauges_in(k) == gauges_manual(:, 2);
    
    gauges_ind(k) = gauges_manual(tmp); %#ok
end

num_gauges     = length(gauges_ind);
gauges_loc     = [LON(gauges_ind); LAT(gauges_ind)] ;
gauges_loc_all = [LON(gauges_manual(:, 1)); LAT(gauges_manual(:, 1))];

xL = -83.27; %round(min(LON), 2); 
xR = -80;    %round(max(LON), 2);
yB = round(min(LAT), 2); 
yT = round(max(LAT), 2);

STATE_FILE    = 'Exp/N80_ol/all_preassim_mean.nc';

x = double(ncread(STATE_FILE, 'qlink1'));

gY = [ 100, 100, 100 ]/255;

BLUE   = [236,231,242; 166,189,219;  43,140,190]/255;
RED    = [254,232,200; 253,187,132; 227, 74, 51]/255;
GREEN  = [229,245,249; 153,216,201;  44,162, 95]/255;

xm = x(:, 240);
lw = 15 * xm/max(xm) + 0.1;

goi = [2297310, 2300017, 2300100]; 

goi_l(1) = find(goi(1) == gauges_manual(:, 2));
goi_l(2) = find(goi(2) == gauges_manual(:, 2));
goi_l(3) = find(goi(3) == gauges_manual(:, 2));

%% 
figure('pos', [2150, 302, 600, 780]); 

xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude', 'FontSize', 16)

set(gca, 'FontSize', 16, 'XLim', [xL, xR], 'YLim', [yB, yT], ...
     'YTick', linspace(yB, yT, 5), ...
     'YTickLabelRotation', 90); grid on; hold on; box on

start_link = 1; 

mval = xm;
cmap = winter(n_links); 

size_cmap = size(cmap, 1); 

zi = xm;
bi = (zi - min(zi))/(max(zi)-min(zi));

cn = ceil(bi * n_links);
cn = max(cn, 1);

for i = start_link:n_links
    % ith link
    lon_i = LON(i);
    lat_i = LAT(i);
    num_i = num_up_links(i);
    col_i = cmap(cn(i), :);  
    wid_i = lw(i);

    if num_i > 0
        z  = fromIndices(fromIndsStart(i):fromIndsEnd(i));
        x1 = lon_i*ones(1, num_i); y1 = LON(z); p1 = [x1; y1];
        x2 = lat_i*ones(1, num_i); y2 = LAT(z); p2 = [x2; y2];
        line(p1, p2, 'color', col_i, 'LineWidth', wid_i)
    end  
end

% borders('states', 'Color', 'k', 'LineWidth', 1.5); grid on
% To display state borders: you can download the borders routine
% from "the climate data toolbox for MATLAB" at: 
% https://www.mathworks.com/matlabcentral/fileexchange/50390-borders

[~, uni_cols] = unique(cn); 
num_colors    = length(uni_cols);
col_levels    = cmap(cn(uni_cols), :);

[~, r] = sort(mval(uni_cols));

colormap(gca, col_levels(r,:));

hc = colorbar; 

hc_cubes   = min(7, num_colors); 
hc_range   = linspace(0, num_colors, hc_cubes);
hc_locate  = hc_range/num_colors;

hc_range(1)  = 1; 
hc_range     = ceil(hc_range); 
hc_full_cols = sort(round(mval(uni_cols), 5));

ind_tmp = find(hc_full_cols == 0);
bob     = hc_range - ind_tmp;

if size(ind_tmp, 1) > 0
    pos_tmp = find(abs(bob) == min(abs(bob))) ;
    pos_tmp = pos_tmp(end);

    hc_range(pos_tmp) = ind_tmp;
end

hc_labels_tmp = hc_full_cols(hc_range);

hc_labels = round(hc_labels_tmp); 

set(hc, 'YTick', hc_locate, 'YTickLabel', hc_labels, ...
    'Position', [.8, .6, .03, .2], 'FontSize', 14)

set(get(hc, 'title'), 'String', 'm^3/s', 'FontSize', 14)

% Assimilated gauges
for k = 1:num_gauges_all 

    if k == goi_l(1)
        plot(gauges_loc_all(1, k), gauges_loc_all(2, k), 'd', ...
            'MarkerFaceColor', GREEN(3, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 16);
    elseif k == goi_l(2)
        plot(gauges_loc_all(1, k), gauges_loc_all(2, k), 's', ...
            'MarkerFaceColor', GREEN(3, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 16);
    elseif k == goi_l(3)
        plot(gauges_loc_all(1, k), gauges_loc_all(2, k), '>', ...
            'MarkerFaceColor', GREEN(3, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 16);
    else
        plot(gauges_loc_all(1, k), gauges_loc_all(2, k), 'o', ...
            'MarkerFaceColor', RED(3, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
end

tt = text(-80.6, 29.07, 'Streamflow', 'FontSize', 20);
set(tt, 'Rotation', 90);

title('WRF-Hydro Flooding Domain', 'FontSize', 20)
