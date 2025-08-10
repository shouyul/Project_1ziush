function [lowest_err, align_sample, align_time, num_valid_samples] = findStreamAlignmentFromCNT(EEG_cnt,...
    data2match, chans, segment_lats, searchWindow, searchDirection, do_plot)
% Inputs:
% EEG_cnt           - data extracted from .cnt file with pop_loadeep_v4_custom
% data2match        - pattern of data to find in the .cnt recording
%                       [chans x samples] matrix
% chans             - over which channels this pattern should be found.
%                       Should match the number of rows in data2match
%                       [vector] of channel indexes, ex: 1:63 for the first
%                       63 channels.
% segment_lats      - [int, int] Latencies in the .cnt recording between which you expect to
%                       find the data2match
% searchWindow      - [min max] time window between which you expect the data2match to
%                       begin (in seconds from the beginning of the .cnt recording)
% searchDirection   - 'forward' or 'backward'
% do_plot           - boolean to indicate if you want to plot the result of
%                       the search.

%% Checks
if length(chans) ~= size(data2match,1)
    error('Number of channnels inconsistent')
end

%% Main
num_samples2match = size(data2match,2);
comps = 0;

if EEG_cnt.srate == 1000
    target = round(sum(searchWindow)*1000/2);
elseif EEG_cnt.srate == 500
    if mod(round(sum(searchWindow)*1000/2),2) == 0
        target = round(sum(searchWindow)*1000/2);
    else
        target = round(sum(searchWindow)*1000/2)-1;
    end
else
    error('To implement: rounding may pose issues depending on your sampling rate')
end

for t = 1:(EEG_cnt.pnts-num_samples2match)
    if EEG_cnt.times(t) - EEG_cnt.times(segment_lats(1)) == target
        ref_t = t;
    end
    if (EEG_cnt.times(t) - EEG_cnt.times(segment_lats(1)) >= searchWindow(1)*1000)...
            && (EEG_cnt.times(t) - EEG_cnt.times(segment_lats(1)) < searchWindow(2)*1000)
        comps = comps+1;
        switch searchDirection
            case 'forward'
                comp_results(comps,1) = t+num_samples2match-1;
            case 'backward'
                comp_results(comps,1) = t;
        end
        % Difference between data
        diff_samples = abs(EEG_cnt.data(chans,t:t+num_samples2match-1) - data2match);
        comp_results(comps,2) = mean(diff_samples,'all');
        comp_results(comps,3) = std(diff_samples,1,'all');
    end
end

[lowest_err, i] = min(comp_results(:,2));
align_sample = comp_results(i,1);
align_time = (EEG_cnt.times(align_sample) - EEG_cnt.times(ref_t))/1000;
switch searchDirection
    case 'forward'
        if segment_lats(1) > (1+align_sample)
            warning('You may not be fitting the right segment')
        end
        num_valid_samples = segment_lats(2) - (1+align_sample);
    case 'backward'
        if segment_lats(2) < (1+align_sample)
            warning('You may not be fitting the right segment')
        end
        num_valid_samples = align_sample - (1 + segment_lats(1));
end

%% Summary plot of the search
if do_plot
    X = EEG_cnt.times(comp_results(:,1))./1000;
    Y = comp_results(:,2)'./1000;
    Y_pos_std = Y + (comp_results(:,3)'./2000);
    Y_neg_std = Y - (comp_results(:,3)'./2000);
    figure;
    hold on
    p = plot(X, Y);
    patch('XData', [X, fliplr(X)], 'YData', [Y_pos_std, fliplr(Y_neg_std)],...
        'FaceColor', p.Color, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    scatter(X(i), lowest_err/1000, 20, 'r', 'filled')
    ylim([0,inf])
    xlabel('Time (s)')
    ylabel('Mean absolute difference (mV)')
    legend({'Mean', 'STD', sprintf('Lowest point: %.3f +/- %.3f mV',...
        lowest_err/1000, comp_results(i,3)/1000)})
end
end

