function plot_distribution(distrib, rejected, datatype, pipeline, threshold, location)
% Plot distributions to enlight the outliers selected by the APP pipeline
% Does not assign xlabel and title as it will vary too often
% Inputs:
% distrib       - vector of the data to plot (histogram)
% outliers      - outliers indices
% datatype      - datatype for the distribution ('Channel' or 'Epoch')
%                   Only useful for y label
% pipeline      - pipeline used to reject the bad datatype ('APP_both', 'APP_left', 'APP_right', 'PREP')
% threshold     - (optional) in case a threshold was used and you want it
%                   on the plot
% location      - (optional) legend location ('northwest', 'northeast')

if ~exist('location','var') %default
    location = 'northeast';
end

distrib = reshape(distrib,length(distrib),1);

%{
span = max(distrib)- min(distrib);
span_order = find_order(span);
min_edge = floor(min(distrib)/power(10,span_order))*power(10,span_order);
max_edge = ceil(max(distrib)/power(10,span_order))*power(10,span_order);
step_edge = power(10,span_order-1);
edges = [min_edge:step_edge:max_edge]';
%}

% Just to get information about the bin counting
figure;
h= histogram(distrib, 'Normalization', 'count');
edges = h.BinEdges';
step_edge = h.BinWidth;
max_count = max(h.Values);
max_order = find_order(max_count);
init_step = max(power(10, max_order-1), 1); % in case order gets 0 (then repmat cannot work)
non_zero_bins = find(h.Values > 0);
close(gcf)
clear h

outlier_bins = find_outlier_bins(edges, distrib(rejected));
distrib_augment = [distrib;repmat(edges(non_zero_bins)+step_edge/2,init_step,1)];
distrib_augment_outliers = [distrib(rejected);repmat(edges(outlier_bins)+step_edge/2,init_step,1)];

figure;
hold on
histogram(distrib_augment, edges, 'Normalization', 'count');
histogram(distrib_augment_outliers, edges, 'Normalization', 'count', 'FaceAlpha', 1);
if exist('threshold','var') && ~isempty(threshold)
    line([threshold threshold], [init_step, max_count+init_step], 'LineStyle', '--', 'Color', 'r')
end
line([(edges(1)+edges(2))/2,(edges(end-1)+edges(end))/2], [init_step, init_step], 'LineStyle', '--', 'Color', 'k')
switch pipeline
    % No threshold in APP
    case 'APP_both'
        legend({'Original distribution','Outliers', 'True zero'}, 'Location', location)
    case 'APP_right'
        legend({'Original distribution','Outliers (right tail only)', 'True zero'}, 'Location', location)
    case 'APP_left'
        legend({'Original distribution','Outliers (left tail only)', 'True zero'}, 'Location', location)
    case 'PREP'
        if exist('threshold','var') && ~isempty(threshold)
            legend({'Original distribution','Rejected', 'Rejection threshold','True zero'}, 'Location', location)
        else
            legend({'Original distribution','Rejected', 'True zero'}, 'Location', location)
        end
    otherwise
        warning('No such pipeline')
end
% Realign the y labels
tick_step = power(10, max_order-1)*5;
ticks2keep = init_step:tick_step:max_count+init_step+tick_step;
yticks(ticks2keep)
ylabels = cell(length(ticks2keep),1);
for t = 1:numel(ylabels)
    ylabels{t} = num2str(ticks2keep(t)-init_step);
end
yticklabels(ylabels)
ylabel([datatype ' count'])

    function order = find_order(value)
        order = 0;
        while value/power(10,order)>10
            order = order +1;
        end
        while value/power(10,order)<1
            order = order-1;
        end
    end

    function bins = find_outlier_bins(edges, outlier_values)
        bins = [];
        for b = 1:length(edges)-1
            if ~isempty(intersect(find(outlier_values >= edges(b)), find(outlier_values < edges(b+1))))
                bins = [bins, b];
            end
        end
    end
end

