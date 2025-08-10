function EEG_app = rejectBadTempsWithAPP(EEG, cfg, plot)
subject = cfg.subjects(cfg.current_subject).id;

% 1. HP filter
lowcutoff = cfg.filterPreProc.low_cut_off;
highcutoff = cfg.filterPreProc.high_cut_off;
if ~isempty(lowcutoff)
    fprintf('Highpass Filtering (%.1f Hz) for APP...\n', lowcutoff)
end
if ~isempty(highcutoff)
    fprintf('Lowpass Filtering (%.1f Hz) for APP...\n', highcutoff)
end
[EEG_HP] = custom_filter(EEG, lowcutoff, highcutoff);

% Add path to prepPipeline subdirectories if not in the list
tmp = which('getPipelineDefaults');
if isempty(tmp)
    myPath = fileparts(which('prepPipeline'));
    addpath(genpath(myPath));
end

% 2. Remove Line Noise with PREP pipeline functions
disp('Removing Line Noise...')
[EEG_HP_noLN, lineNoiseOut] = removeLineNoise_custom(EEG_HP, cfg.lineNoiseRemoval_method, false);
clear EEG_HP
% If you want to save the filteredEEG_noLN struct with the LineNoiseRemoval information for later:
EEG_HP_noLN.etc.lineNoiseRemoval = lineNoiseOut;

% 3. Detect transient artifacts with APP
% Make figures folder if not created yet
if ~exist(fullfile(cfg.figures_folder, 'APP_distributions'), 'dir')
    mkdir(fullfile(cfg.figures_folder, 'APP_distributions'));
end

% Remove all events (simpler for later)
% Additionally, avoid loosing the latency information because of boundaries
% (regepochs will delete any epochs overlapping between boundaries
bound_events = EEG_HP_noLN.event(strcmp({EEG_HP_noLN.event.type}, 'boundary'));
EEG_noevent = pop_editeventvals(EEG_HP_noLN, 'delete', 1:length(EEG_HP_noLN.event));

% Epoch the EEG
disp('APP: Epoching data for Bad segments removal...')
duration = EEG_noevent.xmax; %duration in s
Nmax_wdws = floor(duration/1); % 1 indicates that we aim for 1s epochs
epoch_wdw = duration/Nmax_wdws; % in seconds

EEG_epoched = eeg_regepochs(EEG_noevent, 'recurrence', epoch_wdw, 'limits', [0 epoch_wdw],...
    'rmbase', NaN, 'extractepochs', 'on');
% Some additional events have appeared...
events2del = find([EEG_epoched.event.latency]~=round([EEG_epoched.event.latency]));
EEG_epoched = pop_editeventvals(EEG_epoched,'delete', events2del);
clear EEG_noevent

% Find epochs concerned by discontinuity
if ~isempty(bound_events)
    discountinuous_epochs = [];
    count = 1;
    for b = 1:length(bound_events)
        ep = find([EEG_epoched.event.latency] > bound_events(b).latency, 1);
        if ~isempty(ep) && ep > 1
            discountinuous_epochs(count) = ep-1;
            count = count+1;
        end
    end
end
discountinuous_epochs = unique(discountinuous_epochs);
n_epochs_disc = length(discountinuous_epochs);

% Avoid using the discountinuous epochs in the calculation
epochs2use = setdiff(1:EEG_epoched.trials, discountinuous_epochs);
% Remove EOG channels if any
EOG_channels = strcmp({EEG.chanlocs.type},'EOG');
Data = EEG_epoched.data(~EOG_channels,:,epochs2use);

% Max amplitude difference
Max_amps = squeeze(max(Data, [], 2));
Min_amps = squeeze(min(Data, [], 2));
Amp_diffs = Max_amps - Min_amps;

[Epoch_mean_amp_diff,~] = Biweight_custom(Amp_diffs', EEG_epoched.nbchan - sum(EOG_channels), cfg.APP.censorBiweight);
disp('APP: Detecting bad epochs by Mean amplitude difference criterion')
RejEpochs4AmpDiff = FindOutliers_custom(Epoch_mean_amp_diff, cfg.APP, 'right');
plot_distribution(Epoch_mean_amp_diff, RejEpochs4AmpDiff, 'Epoch', 'APP_right')
xlabel('Mean amplitude difference [microV]')
title({sprintf('%s - Distribution of the mean amplitude difference', subject),...
    'for the epochs inspected by APP'})
saveCurrentFig([cfg.figures_folder 'APP_distributions' filesep],...
    sprintf('%s_epochs_mean_amp_diff', subject), {'png'}, [600 500]);

% Epoch variance or the mean GFP
Epoch_GFP = mean(squeeze(std(Data,0,1)));
disp('APP: Detecting bad epochs by GFP criterion')
RejEpochs4GFP = FindOutliers_custom(Epoch_GFP, cfg.APP, 'right');
plot_distribution(Epoch_GFP, RejEpochs4GFP, 'Epoch', 'APP_right')
xlabel('Global field potential [microV]')
title({sprintf('%s - Distribution of the global field potential', subject),...
    'for the epochs inspected by APP'})
saveCurrentFig([cfg.figures_folder 'APP_distributions' filesep],...
    sprintf('%s_epochs_GFP', subject), {'png'}, [600 500]);

% Epoch's mean deviation from channel means.
sz_2 = (EEG_epoched.trials-n_epochs_disc)*EEG_epoched.pnts;
[Chans_mean,~] = Biweight_custom(reshape(Data, EEG_epoched.nbchan - sum(EOG_channels), sz_2), sz_2,...
    cfg.APP.censorBiweight); % channel mean for all epochs
Epoch_mean_dev = mean(abs(squeeze(mean(Data,2)) - Chans_mean));
disp('APP: Detecting bad epochs by deviation criterion')
RejEpochs4MeanDev = FindOutliers_custom(Epoch_mean_dev, cfg.APP, 'right');
plot_distribution(Epoch_mean_dev, RejEpochs4MeanDev, 'Epoch', 'APP_right')
xlabel('Mean deviation from channel mean [microV]')
title({sprintf('%s - Distribution of the mean deviation', subject),...
    'for the epochs inspected by APP'})
saveCurrentFig([cfg.figures_folder 'APP_distributions' filesep],...
    sprintf('%s_epochs_mean_dev', subject), {'png'}, [600 500]);

% Join the rejected epochs
RejEpochs = sort(union(union(RejEpochs4AmpDiff, RejEpochs4GFP), RejEpochs4MeanDev));
% Add neighbouring epochs to the rejection
RejEpochs = sort(union(union(RejEpochs, RejEpochs-1), RejEpochs+1));
RejEpochs = RejEpochs(RejEpochs>0);
RejEpochs = RejEpochs(RejEpochs<=length(epochs2use));
% Add discontinuous epochs
BadEpochs = sort(union(epochs2use(RejEpochs), discountinuous_epochs));

% Define the rejection intervals:
rejected_intervals = [];
count = 1;
for ep = 1:length(BadEpochs)
    if ep == 1
        rejected_intervals(count,1) = EEG_epoched.event(BadEpochs(ep)).latency;
        rejected_intervals(count,2) = EEG_epoched.event(BadEpochs(ep)).latency + EEG_epoched.pnts - 1;
    elseif BadEpochs(ep-1) == BadEpochs(ep)-1 % contiguous epochs
        % just move the end of the interval
        rejected_intervals(count,2) = EEG_epoched.event(BadEpochs(ep)).latency + EEG_epoched.pnts - 1;
    else
        % new interval
        count = count+1;
        rejected_intervals(count,1) = EEG_epoched.event(BadEpochs(ep)).latency;
        rejected_intervals(count,2) = EEG_epoched.event(BadEpochs(ep)).latency + EEG_epoched.pnts - 1;
    end
end

%         % Define the segment indices to reject:
%         rejected_segments_index = zeros(length(BadEpochs), 2);
%         for i = 1:length(RejEpochs)
%             epoch_local_ind = find([EEG_epoched.event.epoch]==RejEpochs(i));
%             epoch_global_ind = EEG_epoched.event(epoch_local_ind).urevent;
%             rejected_segments_index(i,1) = EEG_epoched.urevent(epoch_global_ind).latency;
%             rejected_segments_index(i,2) = rejected_segments_index(i,1) + EEG_epoched.pnts;
%         end

EEG_app = pop_select(EEG_HP_noLN, 'nopoint', rejected_intervals);

clean_sample_mask = ones(1,EEG_HP_noLN.pnts);
for i = 1:size(rejected_intervals,1)
    clean_sample_mask(rejected_intervals(i,1):rejected_intervals(i,2))=0;
end

if plot
    fprintf('%.1f%% of data kept with APP.\n', 100*sum(clean_sample_mask)/length(clean_sample_mask))
    EEG_app.etc.clean_sample_mask = logical(clean_sample_mask);
    vis_artifacts(EEG_app, EEG_HP_noLN);
    EEG_app.etc = rmfield(EEG_app.etc, 'clean_sample_mask');
    close gcf
end

% Structure saving all parameters:
epoch_cleaning = struct('epoch_wdw', epoch_wdw,...
    'censorBiweigth', cfg.APP.censorBiweight,...
    'z_criterion', cfg.APP.z_criterion,...
    'inner_fence', cfg.APP.inner_fence,...
    'skew_side_limit', cfg.APP.skew_side_limit,...
    'opp_side_limit', cfg.APP.opp_side_limit,...
    'Rejected_epochs', RejEpochs,...
    'epochs_ampDiff', Epoch_mean_amp_diff,...
    'Rejected_epochs_ampDiff', RejEpochs4AmpDiff,...
    'epochs_GFP', Epoch_GFP,...
    'Rejected_epochs_GFP', RejEpochs4GFP,...
    'epochs_MeanDev', Epoch_mean_dev,...
    'Rejected_epochs_MeanDev', RejEpochs4MeanDev,...
    'rejected_byAPP_latencies', rejected_intervals);

EEG_app.etc.APP.epoch_cleaning = epoch_cleaning;
EEG_app.etc.APP.rejectedSamples = ~logical(clean_sample_mask);
end
