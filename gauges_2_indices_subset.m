function obs_ind_id = gauges_2_indices_subset(route)

gauges  = strtrim(ncread(route, 'gages')');
n_links = length(gauges);

k = 0;
for i = 1:n_links
    ob_id = gauges(i, :);
    
    if sum(isspace(ob_id)) < 15
        % found a gauge
        k = k + 1;
        obs_ind_id(k, :) = [i, str2double(ob_id)]; %#ok
    end
end