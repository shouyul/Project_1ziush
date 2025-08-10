function  bridgesOut = findBridgedChannels(signal, bridgesIn, do_plot)
% Inspired from findNoisyChannels function from PREP pipeline
% Inspired from get_elecdists function from Scott Burwell
% https://github.com/sjburwell/eeg_commander/blob/master/functions/basefunc/get_elecdists.m


%% Check the incoming parameters
if nargin < 1
    error('findNoisyChannels:NotEnoughArguments', 'requires at least 1 argument');
elseif isstruct(signal) && ~isfield(signal, 'data')
    error('findNoisyChannels:NoDataField', 'requires a structure data field');
elseif size(signal.data, 3) ~= 1
    error('findNoisyChannels:DataNotContinuous', 'data must be a 2D array');
elseif nargin < 2 || ~exist('bridgesIn', 'var') || isempty(bridgesIn)
    bridgesIn = struct();
end

%% Set the defaults and initialize as needed
bridgesOut = bridgesIn;
% highCorrOut = getHighCorrStructure();
% defaults = getPrepDefaults(signal, 'reference');
% [highCorrOut, errors] = checkPrepDefaults(highCorrIn, highCorrOut, defaults);
% if ~isempty(errors)
%     error('findNoisyChannels:BadParameters', ['|' sprintf('%s|', errors{:})]);
% end
%% Fix the channel locations
channelLocations = bridgesOut.chanlocs;
evaluationChannels = sort(bridgesOut.evaluationChannels); % Make sure channels are sorted
evaluationChannels = evaluationChannels(:)';            % Make sure row vector
bridgesOut.evaluationChannels = evaluationChannels;
originalChannels = 1:size(signal.data, 1);

%% Extract the data required
data = signal.data;
originalNumberChannels = size(data, 1);          % Save the original channels
data = double(data(evaluationChannels, :))';      % Remove the unneeded channels
[signalSize, numberChannels] = size(data);
epochFrames = bridgesOut.epochWindowSeconds * signal.srate;
epochWindow = 0:(epochFrames - 1);
epochOffsets = 1:epochFrames:(signalSize-epochFrames);
Nepochs = length(epochOffsets);

%% Compute metrics
channel_elecDist = nan(numberChannels, numberChannels, Nepochs);
channel_correlations = nan(numberChannels, numberChannels, Nepochs);
dataWin = reshape(data(1:epochFrames*Nepochs, :)', numberChannels, epochFrames, Nepochs);

ppm = ParforProgressbar(Nepochs, 'showWorkerProgress', true,...
    'title', 'Computing channels bridge metrics');
parfor k = 1:Nepochs
    dataPortion = squeeze(dataWin(:, :, k))';
    % Electrical distance
    for ch = 1:numberChannels
        diffdataPortion = repmat(dataPortion(:,ch),[1 numberChannels]) - dataPortion;
        channel_elecDist(ch, :, k) = var(diffdataPortion);
    end
    
    % Correlation
    windowCorrelation = corrcoef(dataPortion);
    abs_corr = abs(windowCorrelation - diag(diag(windowCorrelation)));
    channel_correlations(:, :, k)  = abs_corr;
    %channelCorrelations(k, :)  = quantile(abs_corr, 0.98);
    ppm.increment();
end
delete(ppm);
clear dataWin;

%% Analyse electrical distance
% Following algorithm described in https://www.sciencedirect.com/science/article/pii/S1388245713010158
% Take triangular superior part of ED matrix
for e = 1:Nepochs
    channel_elecDist(:,:,e) = triu(channel_elecDist(:,:,e));
end
channel_elecDist(channel_elecDist == 0) = NaN;
% Scale ED matrix
channel_elecDist = channel_elecDist*(100/median(channel_elecDist(:), 'omitnan'));
% Make Histogram distribution
binSizeHist = 0.25;
edgesHist = 0:binSizeHist:max(channel_elecDist(:));
[counts, ~] = histcounts(channel_elecDist(:),edgesHist);
binSizeInterp = 0.05;
evalHist = edgesHist(1:end-1) + binSizeHist/2;
evalInterp = evalHist(1):binSizeInterp:evalHist(end);
counts = interp1(evalHist,counts,evalInterp,'pchip');
maxED = 10;
if do_plot
    figure
    plot(evalInterp, counts)
    xlim([0,maxED])
end

% Find local miminum with ED<maxED
x_cands = evalInterp < maxED;
[~,xmin] = min(counts(x_cands));
if xmin == 1
    disp('No evidence for electrodes bridging from electrical distance analysis');
    ED_cutoff = NaN;
    lowED_chans = [];
    clusters_ED = {};
else
    % Find bridged channels (if exist)
    ED_cutoff = evalInterp(xmin);
    EDtoolow_allepochs = channel_elecDist <= ED_cutoff;
    EDtoolow_mask = false(numberChannels,numberChannels);
    for ch1 = 1:(numberChannels-1)
        for ch2 = (ch1+1):numberChannels
            if sum(EDtoolow_allepochs(ch1,ch2,:)) > (Nepochs/2)
                fprintf('%s and %s seem bridged\n',...
                    channelLocations(evaluationChannels(ch1)).labels,...
                    channelLocations(evaluationChannels(ch2)).labels)
                EDtoolow_mask(ch1,ch2) = true;
            end
        end
    end
    
    if any(EDtoolow_mask, 'all')
        % Find clusters
        [lowED_chans, clusters_ED] = findClustersInMask(EDtoolow_mask);
        % Remap to original channels
        lowED_chans = evaluationChannels(lowED_chans);
        for cl = 1:numel(clusters_ED)
            clusters_ED{cl} = evaluationChannels(clusters_ED{cl});
        end
    else
        lowED_chans = [];
        clusters_ED = {};
    end
end

bridgesOut.elecDistCutoff = ED_cutoff;
bridgesOut.medianElecDist = median(channel_elecDist,3);
bridgesOut.lowElecDistChans = lowED_chans;
bridgesOut.lowElecDistClusters = clusters_ED;

clear channel_elecDist counts evalHist evalInterp;

%% Analyse correlations
% Take triangular superior part
for e = 1:Nepochs
    channel_correlations(:,:,e) = triu(channel_correlations(:,:,e));
end
channel_correlations(channel_correlations == 0) = NaN;
% Make Histogram distribution
binSizeHist = 0.001;
edgesHist = 0:binSizeHist:1;
[counts, ~] = histcounts(channel_correlations(:),edgesHist);
binSizeInterp = 0.0002;
evalHist = edgesHist(1:end-1) + binSizeHist/2;
evalInterp = 0:binSizeInterp:1;
counts = interp1(evalHist,counts,evalInterp,'pchip');
minCorr = 0.95;
if do_plot
    figure
    plot(evalInterp, counts)
    xlim([minCorr,1])
end
% Find local miminum with corr>minCorr
x_cands = find(evalInterp > minCorr);
[~,xmin] = min(counts(x_cands));
if xmin == length(x_cands)
    disp('No evidence for electrodes bridging from correlation analysis');
    corr_cutoff = NaN;
    lowHC_chans = [];
    clusters_HC = {};
else
    % Find bridged channels (if exist)
    corr_cutoff = evalInterp(x_cands(xmin));
    corrTooHigh_allepochs = channel_correlations >= corr_cutoff;
    corrTooHigh_mask = false(numberChannels,numberChannels);
    for ch1 = 1:(numberChannels-1)
        for ch2 = (ch1+1):numberChannels
            if sum(corrTooHigh_allepochs(ch1,ch2,:)) > (Nepochs/2)
                fprintf('%s and %s seem bridged\n',...
                    channelLocations(evaluationChannels(ch1)).labels,...
                    channelLocations(evaluationChannels(ch2)).labels)
                corrTooHigh_mask(ch1,ch2) = true;
            end
        end
    end
    
    if any(corrTooHigh_mask, 'all')
        % Find clusters
        [lowHC_chans, clusters_HC] = findClustersInMask(corrTooHigh_mask);
        % Remap to original channels
        lowHC_chans = evaluationChannels(lowHC_chans);
        for cl = 1:numel(clusters_HC)
            clusters_HC{cl} = evaluationChannels(clusters_HC{cl});
        end
    else
        lowHC_chans = [];
        clusters_HC = {};
    end
end

% highPercCorr = nan(numberChannels, numberChannels);
% for ch1 = 1:(numberChannels-1)
%     for ch2 = (ch1+1):numberChannels
%         highPercCorr(ch1,ch2) = quantile(squeeze(channel_correlations(ch1,ch2,:)), 1 - bridgesIn.badTimeThreshold);
%         %corr_quant(ch2,ch1) = corr_quant(ch1,ch2);
%     end
% end
%
% highCorr_mask = highPercCorr > bridgesIn.correlationThreshold;
%
% %% Find clusters
% [highCorr_chans, clusters] = findClustersInMask(highCorr_mask);
% % Remap to original channels
% highCorr_chans = evaluationChannels(highCorr_chans);
% for cl = 1:numel(clusters)
%     clusters{cl} = evaluationChannels(clusters{cl});
% end

% channelsBycluster = zeros(1,length(originalChannels));
% for cl = 1:numel(clusters)
%     channelsBycluster(clusters{cl}) = cl;
% end
%
% figure
% topoplot(channelsBycluster*10, channelLocations,...
%     'plotgrid', highCorrIn.gridPattern.grid_inds, 'maplimits', 'absmax');
% title({'Clusters of highly correlated channels',...
%     sprintf('Correlation above %.3f in at least %d%% of the recording',highCorrIn.correlationThreshold,highCorrIn.badTimeThreshold*100),...
%     sprintf('Found %d clusters',numel(clusters)),...
%     'Green = not correlated'});


bridgesOut.corrCutoff = corr_cutoff;
bridgesOut.medianCorr = median(channel_correlations,3);
bridgesOut.highCorrChans = lowHC_chans;
bridgesOut.highCorrClusters = clusters_HC;

% bridgesOut.corrQuantile = highPercCorr;
% bridgesOut.highCorrChans = highCorr_chans;
% bridgesOut.highCorrClusters = clusters;

%% Helper functions for findNoisyChannels
function noisyOut = getHighCorrStructure()
noisyOut = struct('srate', [], ...
    'samples', [], ...
    'evaluationChannels', [], ...
    'channelLocations', [], ...
    'correlationWindowSeconds', [], ...
    'correlationThreshold', [], ...
    'badTimeThreshold', []);

