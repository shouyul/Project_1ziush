clear all;
config_PIONEER_Tumbler;
subject_inds = [1];

eyeDetectionChannels = {'R3Z'};
%eyeDetectionChannels = {'VEOGR','Z2Z','Z3Z','R3Z','L3Z'};
plot_kurtosis = false;

for subject_ind = subject_inds
    % Launch eeglab
    if ~exist('EEG', 'var')
        launchEEGLAB
    end
    
    % Overwrite subject for testing
    %subject_ind = 1;
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    EEG = pop_loadset('filename', N.preparedFile,'filepath', N.searchFolder_2);
    %EEG = pop_loadset('filename', N.IClabelledFile,'filepath', N.searchFolder_2arch_rej);
    %EEG.etc = rmfield(EEG.etc, 'filter');
    
    %% Eyes opening
    EEG_sel = pop_select(EEG, 'channel', eyeDetectionChannels);
    EEG_filt1 = custom_filter(EEG_sel, 0.75, 5);
    EEG_epoch1 = pop_epoch(EEG_filt1, {'OpenEyesSoundOn'}, [-10,26]);
    
    EEG_filt2 = custom_filter(EEG_sel, 0.05, 5);
    EEG_epoch2 = pop_epoch(EEG_filt2, {'OpenEyesSoundOn'}, [-10,26]);
    
    % Parameters
    quant = 0.9;
    percMax4Prom = 0;
    wdw = 0.4; % in seconds, for diff calculation
    wdw_tol_bef = 1; % in seconds, for latency criteria
    wdw_tol_aft = 5; % in seconds, for latency criteria
    diffminEO = 5;
    diffmaxEC = -10;
    
    refLatEO = find(EEG_epoch2.times == 0);
    refLatEC = find(EEG_epoch2.times == 20860);
    
    varTypes = {'double','double','double','double','double','double','double','double','double','double'};
    varNames = {'Max', sprintf('Quantile%d',quant*100),'Npeaks','Rank','Latency',...
        'Height1','Width1','Prominence1','Height2','Diff2'};
    EyeOpening = table('Size',[EEG_epoch1.trials,numel(varNames)], 'VariableTypes',varTypes, 'VariableNames', varNames);
    EyeClosing = table('Size',[EEG_epoch1.trials,numel(varNames)], 'VariableTypes',varTypes, 'VariableNames', varNames);
    
    for tr = 1:EEG_epoch1.trials
        EEGmax = max(EEG_epoch1.data(:,:,tr), [], 'all');
        EEGperc = quantile(EEG_epoch1.data(:,:,tr), quant, 'all');
        [heights, locs, widths, proms] = findpeaks(squeeze(EEG_epoch1.data(:,:,tr)),...
            'MinPeakHeight', EEGperc, 'MinPeakProminence', EEGmax*percMax4Prom);
        
        toExclude = locs < wdw*EEG.srate | locs > EEG_epoch1.pnts - wdw*EEG.srate;
        heights = heights(~toExclude);
        locs = locs(~toExclude);
        widths = widths(~toExclude);
        proms = proms(~toExclude);
        
        Npeaks = length(locs);
        EyeOpening.Max(tr) = EEGmax;
        EyeOpening.(sprintf('Quantile%d',quant*100))(tr) = EEGperc;
        EyeOpening.Npeaks(tr) = Npeaks;
        EyeClosing.Max(tr) = EEGmax;
        EyeClosing.(sprintf('Quantile%d',quant*100))(tr) = EEGperc;
        EyeClosing.Npeaks(tr) = Npeaks;
        
        % Default values
        EyeOpening.Rank(tr) = NaN;
        EyeOpening.Latency(tr) = 1;
        EyeOpening.Height1(tr) = NaN;        
        EyeOpening.Width1(tr) = NaN;
        EyeOpening.Prominence1(tr) = NaN;
        EyeOpening.Height2(tr) = NaN;  
        EyeOpening.Diff2(tr) = NaN;
        
        EyeClosing.Rank(tr) = NaN;
        EyeClosing.Latency(tr) = 1;
        EyeClosing.Height1(tr) = NaN;        
        EyeClosing.Width1(tr) = NaN;
        EyeClosing.Prominence1(tr) = NaN;
        EyeClosing.Height2(tr) = NaN;  
        EyeClosing.Diff2(tr) = NaN;
        
        if Npeaks > 1
            latsLeft = -round(wdw*EEG.srate):-1;
            meanLeft = nan(1,Npeaks);
            latsRight = 1:round(wdw*EEG.srate);
            meanRight = nan(1,Npeaks);
            for p = 1:Npeaks
                meanLeft(p) = mean(EEG_epoch2.data(:,locs(p)+latsLeft,tr),'all');
                meanRight(p) = mean(EEG_epoch2.data(:,locs(p)+latsRight,tr),'all');
            end
            diffLR = meanLeft-meanRight;
            
            [diffLR, sortorder] = sort(diffLR, 'descend');
            for p = 1:Npeaks
                if diffLR(p) > diffminEO
                    if locs(sortorder(p)) > (refLatEO - wdw_tol_bef*EEG.srate) &&...
                            locs(sortorder(p)) < (refLatEO + wdw_tol_aft*EEG.srate)
                        EyeOpening.Rank(tr) = p;
                        EyeOpening.Latency(tr) = locs(sortorder(p));
                        EyeOpening.Height1(tr) = heights(sortorder(p));                        
                        EyeOpening.Width1(tr) = widths(sortorder(p));
                        EyeOpening.Prominence1(tr) = proms(sortorder(p));
                        EyeOpening.Height2(tr) = mean(EEG_epoch2.data(:,locs(sortorder(p)),tr),'all');
                        EyeOpening.Diff2(tr) = diffLR(p);
                        break
                    end
                end
            end
            
            %             if tr==7
            %                 0;
            %             end
            for p = Npeaks:-1:1
                if diffLR(p) < diffmaxEC
                    if locs(sortorder(p)) > (refLatEC - wdw_tol_bef*EEG.srate) &&...
                            locs(sortorder(p)) < (refLatEC + wdw_tol_aft*EEG.srate)
                        EyeClosing.Rank(tr) = 1 + Npeaks - p;
                        EyeClosing.Latency(tr) = locs(sortorder(p));
                        EyeClosing.Height1(tr) = heights(sortorder(p));                        
                        EyeClosing.Width1(tr) = widths(sortorder(p));
                        EyeClosing.Prominence1(tr) = proms(sortorder(p));
                        EyeClosing.Height2(tr) = mean(EEG_epoch2.data(:,locs(sortorder(p)),tr),'all');
                        EyeClosing.Diff2(tr) = diffLR(p);
                        break
                    end
                end
            end
            
            %             [val,bestEO] = max(diffLR);
            %             if val > diffminEO
            %                 EyeOpening.Height(tr) = heights(bestEO);
            %                 EyeOpening.Latency(tr) = locs(bestEO);
            %                 EyeOpening.Width(tr) = widths(bestEO);
            %                 EyeOpening.Prominence(tr) = proms(bestEO);
            %                 EyeOpening.Diff(tr) = diffLR(bestEO);
            %             end
            %
            %             [val,bestEC] = min(diffLR);
            %             if val < diffminEC
            %                 EyeClosing.Height(tr) = heights(bestEC);
            %                 EyeClosing.Latency(tr) = locs(bestEC);
            %                 EyeClosing.Width(tr) = widths(bestEC);
            %                 EyeClosing.Prominence(tr) = proms(bestEC);
            %                 EyeClosing.Diff(tr) = diffLR(bestEC);
            %             end
        elseif Npeaks == 1
            meanLeft = mean(EEG_epoch2.data(:,locs+latsLeft,tr),'all');
            meanRight = mean(EEG_epoch2.data(:,locs+latsRight,tr),'all');
            diffLR = meanLeft-meanRight;
            if diffLR > diffminEO
                EyeOpening.Latency(tr) = locs;
                EyeOpening.Height1(tr) = heights;                
                EyeOpening.Width1(tr) = widths;
                EyeOpening.Prominence1(tr) = proms;
                EyeOpening.Height2(tr) = mean(EEG_epoch2.data(:,locs,tr),'all');
                EyeOpening.Diff2(tr) = diffLR;
            elseif diffLR < diffmaxEC
                EyeClosing.Latency(tr) = locs;
                EyeClosing.Height1(tr) = heights;                
                EyeClosing.Width1(tr) = widths;
                EyeClosing.Prominence1(tr) = proms;
                EyeClosing.Height2(tr) = mean(EEG_epoch2.data(:,locs,tr),'all');
                EyeClosing.Diff2(tr) = diffLR;
            end
        end
    end
    
    %% Some stats
    rejTrials = false(size(EyeOpening,1),1);
    % Outliers in max value are probably artifacts
    params.censorBiweight = 6;
    params.z_criterion = 5;
    params.inner_fence = 1.5;
    params.skew_side_limit = 3;
    params.opp_side_limit = 4;
    %     artifactTrials = FindOutliers_custom(bestEyeOpening.Max,params,'right');
    %     rejTrials(artifactTrials) = true;
    
    rejTrialsEO = rejTrials;
    % reject trials when no value was found
    rejTrialsEO(isnan(EyeOpening.Rank)) = true;
    rejTrialsEO(FindOutliers_custom(EyeOpening.Height2,params,'right')) = true;
    %rejTrialsEO(FindOutliers_custom(EyeOpening.Prominence,params,'right')) = true;
    %     % reject trials when too far from the expected latency
    %     rejTrialsEO(EyeOpening.Latency < (refLatEO - wdw_tol_bef*EEG.srate) |...
    %         EyeOpening.Latency > (refLatEO + wdw_tol_aft*EEG.srate)) = true;
    EyeOpening = addvars(EyeOpening,rejTrialsEO,'NewVariableNames',{'Rejected'});
    
    figure;
    subplot(5,2,1)
    hold on;
    edges = 0:5:max(EyeOpening.Height1)+5;
    histogram(EyeOpening.Height1(~EyeOpening.Rejected),edges)
    histogram(EyeOpening.Height1(EyeOpening.Rejected),edges)
    title({'Best Eye Opening', 'Height1'})
    subplot(5,2,3)
    hold on;
    edges = 0:10:max(EyeOpening.Width1)+10;
    histogram(EyeOpening.Width1(~EyeOpening.Rejected),edges)
    histogram(EyeOpening.Width1(EyeOpening.Rejected),edges)
    title('Width1')
    subplot(5,2,5)
    hold on;
    edges = 0:5:max(EyeOpening.Prominence1)+5;
    histogram(EyeOpening.Prominence1(~EyeOpening.Rejected),edges)
    histogram(EyeOpening.Prominence1(EyeOpening.Rejected),edges)
    title('Prominence1')
    subplot(5,2,7)
    hold on;
    edges = 0:10:max(EyeOpening.Height2)+10;
    histogram(EyeOpening.Height2(~EyeOpening.Rejected),edges)
    histogram(EyeOpening.Height2(EyeOpening.Rejected),edges)
    title('Height2') 
    subplot(5,2,9)    
    hold on;
    edges = 0:5:max(EyeOpening.Diff2)+5;
    histogram(EyeOpening.Diff2(~EyeOpening.Rejected),edges)
    histogram(EyeOpening.Diff2(EyeOpening.Rejected),edges)
    title('Difference2')    
    
    rejTrialsEC = rejTrials;
    % reject trials when no value was found
    rejTrialsEC(isnan(EyeClosing.Rank)) = true;
    %rejTrialsEC(FindOutliers_custom(EyeClosing.Height2,params,'right')) = true;
    %rejTrialsEC(FindOutliers_custom(EyeClosing.Prominence,params,'right')) = true;
    %     % reject trials when too far from the expected latency
    %     rejTrialsEC(EyeClosing.Latency < (refLatEC - wdw_tol_bef*EEG.srate) |...
    %         EyeClosing.Latency > (refLatEC + wdw_tol_aft*EEG.srate)) = true;
    EyeClosing = addvars(EyeClosing,rejTrialsEC,'NewVariableNames',{'Rejected'});
    
    subplot(5,2,2)
    hold on;
    edges = 0:5:max(EyeClosing.Height1)+5;
    histogram(EyeClosing.Height1(~EyeClosing.Rejected),edges)
    histogram(EyeClosing.Height1(EyeClosing.Rejected),edges)
    title({'Best Eye Closing', 'Height1'})
    subplot(5,2,4)
    hold on;
    edges = 0:10:max(EyeClosing.Width1)+10;
    histogram(EyeClosing.Width1(~EyeClosing.Rejected),edges)
    histogram(EyeClosing.Width1(EyeClosing.Rejected),edges)
    title('Width1')
    subplot(5,2,6)
    hold on;
    edges = 0:5:max(EyeClosing.Prominence1)+5;
    histogram(EyeClosing.Prominence1(~EyeClosing.Rejected),edges)
    histogram(EyeClosing.Prominence1(EyeClosing.Rejected),edges)
    title('Prominence1')
    subplot(5,2,8)
    hold on;
    edges = 0:10:max(EyeClosing.Height2)+10;
    histogram(EyeClosing.Height2(~EyeClosing.Rejected),edges)
    histogram(EyeClosing.Height2(EyeClosing.Rejected),edges)
    title('Height2')
    subplot(5,2,10)
    hold on;
    edges = (min(EyeClosing.Diff2)-10):10:0;
    histogram(EyeClosing.Diff2(~EyeClosing.Rejected),edges)
    histogram(EyeClosing.Diff2(EyeClosing.Rejected),edges)
    title('Difference2')
    suptitle(subject)
    
    %% Plotting
    % Vizualization trial by trial
    % figure;topoplot(1:127, EEG.chanlocs, 'style', 'blank', 'electrodes', 'labels');
    % Multiple channels at the same time
    %     for tr = 1:EEG_epoch.trials
    %         pop_plotdata(EEG_epoch, 1, 1:5, tr, sprintf('Trial %d',tr), 1, 1, [0,0]);
    %     end
    
    for ch = 1:numel(eyeDetectionChannels)
        % All trials for one channel
        step = 0.25; % in mV
        meanData1 = squeeze(mean(EEG_epoch1.data(ch,:,:),1))/1000; % convert to mV
        meanData2 = squeeze(mean(EEG_epoch2.data(ch,:,:),1))/1000; % convert to mV
        decMeanData1 = zeros(size(meanData1));
        decMeanData2 = zeros(size(meanData2));
        if plot_kurtosis
            kurtData = nan(size(meanData2));
            wdw = 1000; %in ms
            
            for s = 1:EEG_epoch2.pnts
                if EEG_epoch2.times(1)<=EEG_epoch2.times(s)-wdw/2 && EEG_epoch2.times(end)>=EEG_epoch2.times(s)+wdw/2
                    lats = EEG_epoch2.times>=EEG_epoch2.times(s)-wdw/2 & EEG_epoch2.times<=EEG_epoch2.times(s)+wdw/2;
                    kurtData(s,:) = kurtosis(meanData2(lats,:),1,1);
                end
            end
            
            decKurtData = zeros(size(meanData2));
        end
        
        eyeOpPeak4plot = nan(EEG_epoch2.trials,1);
        eyeClPeak4plot = nan(EEG_epoch2.trials,1);
        for tr = 1:EEG_epoch2.trials
            decMeanData1(:,tr) = meanData1(:,tr) - step*(tr-1);
            decMeanData2(:,tr) = meanData2(:,tr) - step*(tr-1);
            if EyeOpening.Latency(tr) > 1
                eyeOpPeak4plot(tr) = meanData2(EyeOpening.Latency(tr),tr) - step*(tr-1);
            end
            if EyeClosing.Latency(tr) > 1
                eyeClPeak4plot(tr) = meanData2(EyeClosing.Latency(tr),tr) - step*(tr-1);
            end
            if plot_kurtosis
                decKurtData(:,tr) = kurtData(:,tr)/1000 - step*(tr-1); % rescale to plot on the same scale
            end
        end
        
        figure;hold on;
        plot(EEG_epoch1.times,decMeanData1,':k');
        plot(EEG_epoch2.times,decMeanData2);
        scatter(EEG_epoch2.times(EyeOpening.Latency(~EyeOpening.Rejected)),...
            eyeOpPeak4plot(~EyeOpening.Rejected),'g','o')
        scatter(EEG_epoch2.times(EyeOpening.Latency(EyeOpening.Rejected)),...
            eyeOpPeak4plot(EyeOpening.Rejected),'r','o')
        scatter(EEG_epoch2.times(EyeClosing.Latency(~EyeClosing.Rejected)),...
            eyeClPeak4plot(~EyeClosing.Rejected),'g','x')
        scatter(EEG_epoch2.times(EyeClosing.Latency(EyeClosing.Rejected)),...
            eyeClPeak4plot(EyeClosing.Rejected),'r','x')
        
        if plot_kurtosis
            plot(EEG_epoch2.times,decKurtData,'--');
        end
        xline(0, 'Label', 'Eye Opening event');
        xline(20850, 'Label', 'Eye Closing event');
        for tr = 1:EEG_epoch2.trials
            yline(-step*(tr-1)+EyeOpening.(sprintf('Quantile%d',quant*100))(tr)/1000, ':');
            yl = yline(-step*(tr-1), 'Label',sprintf('Trial%d',tr));
            yl.LabelHorizontalAlignment = 'center';
            yl.LabelVerticalAlignment = 'bottom';
            %yline(-step*(tr-1)-viz_thresh/1000, ':');
        end
        xlabel('Time (ms)');
        ylabel('Amplitude (mV)');
        title(sprintf('%s - %s', subject, EEG_epoch2.chanlocs(ch).labels))
    end
    
    
    
    
    % Blinker toolbox
    %     params = getBlinkerDefaults();
    %     [~, com, blinks, blinkFits, blinkProperties, ...
    %         blinkStatistics, params] = pop_blinker(EEG_prepared);
    %
    %     figure; hold on;
    %     for b = 1:size(blinks.signalData(7).blinkPositions,2)
    %         patch('XData', [blinks.signalData(7).blinkPositions(:,b)',fliplr(blinks.signalData(7).blinkPositions(:,b)')],...
    %             'YData', [1,1,-1,-1],'FaceColor',[0.5,0.5,0.5], 'EdgeColor', 'none');
    %     end
    %     plot(blinks.signalData(7).signal);
end