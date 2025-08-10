clear all;
config_PIONEER_Tumbler;

eyeDetectionChannel = {'VEOGR'};
%eyeDetectionChannels = {'VEOGR','Z2Z','Z3Z','R3Z','L3Z'};
plot_single_level = false;
plot_stats = true;
plot_general_level = true;

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
    
    %% Prepare data (following Sharma et al. 2020 method)
    EEG_sel = pop_select(EEG, 'channel', eyeDetectionChannel);
    if strcmp(eyeDetectionChannel{1},'VEOGR')
        % Invert polarity for EOG channel
        EEG_sel.data = -EEG_sel.data; 
    end
    EEG_filt1 = custom_filter(EEG_sel, 0.75, 5);
    EEG_epoch1 = pop_epoch(EEG_filt1, {'OpenEyesSoundOn'}, [-10,26]);
    
    EEG_filt2 = custom_filter(EEG_sel, 0.05, 5);
    EEG_epoch2 = pop_epoch(EEG_filt2, {'OpenEyesSoundOn'}, [-10,26]);
    
    % General parameters
    artifact_thresh = 500; % in microV
    quant = 97.5;
    percMax4Prom = 0;
    wdw_diff = 0.2; % in seconds, for diff calculation
    diffminEO = 10;
    diffmaxEC = -10;
    latsLeft = -round(wdw_diff*EEG.srate):-1;
    latsRight = 1:round(wdw_diff*EEG.srate);
    
    wdw_lat_bef = 1; % in seconds, for latency criteria
    wdw_lat_aft = 5; % in seconds, for latency criteria
    refLatEO = find(EEG_epoch2.times == 0);
    refLatEC = find(EEG_epoch2.times == 20860);
    
    varNames = {'Trial','Quantile','Latency','Max1','Ave1',...
        'Height1','Width1','Prominence1','Max2','Ave2','Height2','Diff2',...
        'EyeOpening','EyeClosing','Blink','Rejected'};
    varTypes = {'int32','double','int32','double','double',...
        'double','double','double','double','double','double','double',...
        'logical','logical','logical','logical'};
    
    Peaks = table('Size',[0,numel(varNames)], 'VariableTypes',varTypes, 'VariableNames', varNames);
    i = 0;
    
    for tr = 1:EEG_epoch1.trials
        EEGave1 = mean(EEG_epoch1.data(:,:,tr), 'all');
        EEGmax1 = max(EEG_epoch1.data(:,:,tr), [], 'all');
        EEGquant = quantile(EEG_epoch1.data(:,:,tr), quant/100, 'all');
        [heights, locs, widths, proms] = findpeaks(squeeze(EEG_epoch1.data(:,:,tr)),...
            'MinPeakHeight', EEGquant, 'MinPeakProminence', EEGmax1*percMax4Prom);
        
        EEGave2 = mean(EEG_epoch2.data(:,:,tr), 'all');
        EEGmax2 = max(EEG_epoch2.data(:,:,tr), [], 'all');
        
        Npeaks = length(locs);
        for p = 1:Npeaks
            i = i+1;
            Peaks.Trial(i) = tr;
            Peaks.Quantile(i) = quant;
            Peaks.Latency(i) = locs(p);
            Peaks.Max1(i) = EEGmax1;
            Peaks.Ave1(i) = EEGave1;
            Peaks.Height1(i) = heights(p);
            Peaks.Width1(i) = widths(p);
            Peaks.Prominence1(i) = proms(p);
            Peaks.Max2(i) = EEGmax2;
            Peaks.Ave2(i) = EEGave2;
            
            % Default values
            Peaks.EyeOpening(i) = false;
            Peaks.EyeClosing(i) = false;
            Peaks.Blink(i) = false;
            Peaks.Rejected(i) = false;
            if locs(p) > wdw_diff*EEG.srate && locs(p) < EEG_epoch1.pnts - wdw_diff*EEG.srate
                Peaks.Height2(i) = mean(EEG_epoch2.data(:,locs(p),tr),'all');
                
                if Peaks.Height2(i) >= artifact_thresh
                    Peaks.Diff2(i) = NaN;
                    Peaks.Rejected(i) = true;
                    continue
                end
                
                meanLeft = mean(EEG_epoch2.data(:,locs(p)+latsLeft,tr),'all');
                meanRight = mean(EEG_epoch2.data(:,locs(p)+latsRight,tr),'all');
                diffLR = meanLeft-meanRight;
                Peaks.Diff2(i) = diffLR;
                
                if diffLR > diffminEO
                    Peaks.EyeOpening(i) = true;
                elseif diffLR < diffmaxEC
                    Peaks.EyeClosing(i) = true;
                else
                    Peaks.Blink(i) = true;
                end
            else
                Peaks.Height2(i) = NaN;
                Peaks.Diff2(i) = NaN;
                Peaks.Rejected(i) = true;
            end
        end
        
        % Look at the sequencing
        tr_peaks = Peaks.Trial == tr;
        eo_peaks = find(Peaks.EyeOpening & tr_peaks);
        ec_peaks = find(Peaks.EyeClosing & tr_peaks);
        
        Peaks.Rejected(eo_peaks) = true;
        Peaks.Rejected(ec_peaks) = true;
        [diff_eo,sortorder_eo] = sort(Peaks.Diff2(eo_peaks),'descend');
        [diff_ec,sortorder_ec] = sort(Peaks.Diff2(ec_peaks),'ascend');
        
        lats = 1:EEG_epoch1.pnts;
        eo_period = false(1,length(lats));
        c_max = max([length(eo_peaks),length(ec_peaks)]);
        for c = 1:c_max
            if c == 1
                if length(eo_peaks) >= 1 && length(ec_peaks) >= 1
                    Peaks.Rejected(eo_peaks(sortorder_eo(c))) = false;
                    Peaks.Rejected(ec_peaks(sortorder_ec(c))) = false;
                    if Peaks.Latency(eo_peaks(sortorder_eo(c))) < Peaks.Latency(ec_peaks(sortorder_ec(c)))
                        eo_period(lats >= Peaks.Latency(eo_peaks(sortorder_eo(c))) &...
                            lats < Peaks.Latency(ec_peaks(sortorder_ec(c)))) = true;
                    else
                        eo_period(lats >= Peaks.Latency(eo_peaks(sortorder_eo(c))) |...
                            lats < Peaks.Latency(ec_peaks(sortorder_ec(c)))) = true;
                    end
                elseif length(eo_peaks) >= 1
                    Peaks.Rejected(eo_peaks(sortorder_eo(c))) = false;
                    eo_period(lats >= Peaks.Latency(eo_peaks(sortorder_eo(c)))) = true;
                elseif length(ec_peaks) >= 1
                    Peaks.Rejected(ec_peaks(sortorder_ec(c))) = false;
                    eo_period(lats < Peaks.Latency(ec_peaks(sortorder_ec(c)))) = true;
                else
                    eo_period(lats) = true;
                end
            else
                break
                % Not completed yet
                if c <= length(eo_peaks) && ec_period(Peaks.Latency(eo_peaks(sortorder_eo(c))))
                    Peaks.Rejected(eo_peaks(sortorder_eo(c))) = false;
                    
                end
            end
        end
        
        %         [~,best_eo] = max(Peaks.Diff2(eo_peaks));
        %         Peaks.Rejected(eo_peaks(best_eo)) = false;
        %         [~,best_ec] = min(Peaks.Diff2(ec_peaks));
        %         Peaks.Rejected(ec_peaks(best_ec)) = false;
        
        bl_peaks = find(Peaks.Blink & tr_peaks);
        Peaks.Rejected(bl_peaks(~eo_period(Peaks.Latency(bl_peaks)))) = true;
        %
        %         if Peaks.Latency(eo_peaks(best_eo)) < Peaks.Latency(ec_peaks(best_ec))
        %             Peaks.Rejected(bl_peaks(Peaks.Latency(bl_peaks) < Peaks.Latency(eo_peaks(best_eo)) |...
        %                 Peaks.Latency(bl_peaks) > Peaks.Latency(ec_peaks(best_ec)))) = true;
        %         else
        %             Peaks.Rejected(bl_peaks(Peaks.Latency(bl_peaks) < Peaks.Latency(eo_peaks(best_eo)) &...
        %                 Peaks.Latency(bl_peaks) > Peaks.Latency(ec_peaks(best_ec)))) = true;
        %         end
        
        if plot_single_level
            % Plotting
            figure;hold on;
            plot(EEG_epoch1.times,mean(EEG_epoch1.data(:,:,tr),1),':k');
            plot(EEG_epoch2.times,mean(EEG_epoch2.data(:,:,tr),1));
            scatter(EEG_epoch2.times(Peaks.Latency(Peaks.EyeOpening & ~Peaks.Rejected & tr_peaks)),...
                Peaks.Height2(Peaks.EyeOpening & ~Peaks.Rejected & tr_peaks),'g','o')
            scatter(EEG_epoch2.times(Peaks.Latency(Peaks.EyeOpening & Peaks.Rejected & tr_peaks)),...
                Peaks.Height2(Peaks.EyeOpening & Peaks.Rejected & tr_peaks),'r','o')
            scatter(EEG_epoch2.times(Peaks.Latency(Peaks.EyeClosing & ~Peaks.Rejected & tr_peaks)),...
                Peaks.Height2(Peaks.EyeClosing & ~Peaks.Rejected & tr_peaks),'g','x')
            scatter(EEG_epoch2.times(Peaks.Latency(Peaks.EyeClosing & Peaks.Rejected & tr_peaks)),...
                Peaks.Height2(Peaks.EyeClosing & Peaks.Rejected & tr_peaks),'r','x')
            scatter(EEG_epoch2.times(Peaks.Latency(Peaks.Blink & ~Peaks.Rejected & tr_peaks)),...
                Peaks.Height2(Peaks.Blink & ~Peaks.Rejected & tr_peaks),'g','*')
            scatter(EEG_epoch2.times(Peaks.Latency(Peaks.Blink & Peaks.Rejected & tr_peaks)),...
                Peaks.Height2(Peaks.Blink & Peaks.Rejected & tr_peaks),'r','*')
            xline(0, 'Label', 'Eye Opening event');
            xline(20860, 'Label', 'Eye Closing event');
            yline(0,'--');
            yline(EEGquant,':');
            xlabel('Time (ms)');
            ylabel('Amplitude (microV)');
            title(sprintf('Trial %d', tr));
        end
        
    end
    
    if plot_stats
        %% Some stats
        % Outliers in max value are probably artifacts
        params.censorBiweight = 6;
        params.z_criterion = 5;
        params.inner_fence = 1.5;
        params.skew_side_limit = 3;
        params.opp_side_limit = 4;
        %     artifactTrials = FindOutliers_custom(bestEyeOpening.Max,params,'right');
        %     rejTrials(artifactTrials) = true;
        
        figure;
        subplot(3,2,1)
        hold on;
        edges = 0:5:max(Peaks.Height1)+5;
        histogram(Peaks.Height1(Peaks.EyeOpening & ~Peaks.Rejected),edges);
        histogram(Peaks.Height1(Peaks.EyeClosing & ~Peaks.Rejected),edges);
        histogram(Peaks.Height1(Peaks.Blink & ~Peaks.Rejected),edges);
        %histogram(Peaks.Height1(Peaks.Rejected),edges);
        title('Height1')
        subplot(3,2,3)
        hold on;
        edges = 0:10:max(Peaks.Width1)+10;
        histogram(Peaks.Width1(Peaks.EyeOpening & ~Peaks.Rejected),edges);
        histogram(Peaks.Width1(Peaks.EyeClosing & ~Peaks.Rejected),edges);
        histogram(Peaks.Width1(Peaks.Blink & ~Peaks.Rejected),edges);
        %histogram(Peaks.Width1(Peaks.Rejected),edges);
        title('Width1')
        subplot(3,2,5)
        hold on;
        edges = 0:5:max(Peaks.Prominence1)+5;
        histogram(Peaks.Prominence1(Peaks.EyeOpening & ~Peaks.Rejected),edges);
        histogram(Peaks.Prominence1(Peaks.EyeClosing & ~Peaks.Rejected),edges);
        histogram(Peaks.Prominence1(Peaks.Blink & ~Peaks.Rejected),edges);
        %histogram(Peaks.Prominence1(Peaks.Rejected),edges);
        title('Prominence1')
        subplot(3,2,2)
        hold on;
        edges = 0:10:max(Peaks.Height2)+10;
        histogram(Peaks.Height2(Peaks.EyeOpening & ~Peaks.Rejected),edges);
        histogram(Peaks.Height2(Peaks.EyeClosing & ~Peaks.Rejected),edges);
        histogram(Peaks.Height2(Peaks.Blink & ~Peaks.Rejected),edges);
        %histogram(Peaks.Height2(Peaks.Rejected),edges);
        title('Height2')
        subplot(3,2,4)
        hold on;
        edges = (min(Peaks.Diff2)-10):10:max(Peaks.Diff2)+10;
        h1 = histogram(Peaks.Diff2(Peaks.EyeOpening & ~Peaks.Rejected),edges);
        h2 = histogram(Peaks.Diff2(Peaks.EyeClosing & ~Peaks.Rejected),edges);
        h3 = histogram(Peaks.Diff2(Peaks.Blink & ~Peaks.Rejected),edges);
        %h4 = histogram(Peaks.Diff2(Peaks.Rejected),edges);
        title('Difference2')
        %lgd = legend([h1,h2,h3,h4],{'EO','EC','Blink','Rejected'});
        lgd = legend([h1,h2,h3],{'EO','EC','Blink'});
        lgd.Position = [0.6,0.2,0.075,0.1];
        suptitle(subject)
    end
    
    %% Plotting        
    if plot_general_level
        % All trials at once
        step = 0.25; % in mV
        meanData1 = squeeze(mean(EEG_epoch1.data(:,:,:),1))/1000; % convert to mV
        meanData2 = squeeze(mean(EEG_epoch2.data(:,:,:),1))/1000; % convert to mV
        decMeanData1 = zeros(size(meanData1));
        decMeanData2 = zeros(size(meanData2));        

        for tr = 1:EEG_epoch2.trials
            decMeanData1(:,tr) = meanData1(:,tr) - step*(tr-1);
            decMeanData2(:,tr) = meanData2(:,tr) - step*(tr-1);
        end        
        
        peaks4plot = nan(size(Peaks,1),1);
        for p = 1:size(Peaks,1)
            tr = Peaks.Trial(p);
            peaks4plot(p) = decMeanData2(Peaks.Latency(p),tr);
        end
        
        figure;hold on;
        plot(EEG_epoch1.times,decMeanData1,':k');
        plot(EEG_epoch2.times,decMeanData2);
        scatter(EEG_epoch2.times(Peaks.Latency(Peaks.EyeOpening & ~Peaks.Rejected)),...
            peaks4plot(Peaks.EyeOpening & ~Peaks.Rejected),'g','o')
        scatter(EEG_epoch2.times(Peaks.Latency(Peaks.EyeOpening & Peaks.Rejected)),...
            peaks4plot(Peaks.EyeOpening & Peaks.Rejected),'r','o')
        scatter(EEG_epoch2.times(Peaks.Latency(Peaks.EyeClosing & ~Peaks.Rejected)),...
            peaks4plot(Peaks.EyeClosing & ~Peaks.Rejected),'g','x')
        scatter(EEG_epoch2.times(Peaks.Latency(Peaks.EyeClosing & Peaks.Rejected)),...
            peaks4plot(Peaks.EyeClosing & Peaks.Rejected),'r','x')
                scatter(EEG_epoch2.times(Peaks.Latency(Peaks.Blink & ~Peaks.Rejected)),...
            peaks4plot(Peaks.Blink & ~Peaks.Rejected),'g','*')
        scatter(EEG_epoch2.times(Peaks.Latency(Peaks.Blink & Peaks.Rejected)),...
            peaks4plot(Peaks.Blink & Peaks.Rejected),'r','*')
  
        xline(0, 'Label', 'Eye Opening event');
        xline(20850, 'Label', 'Eye Closing event');
        for tr = 1:EEG_epoch2.trials
            yl = yline(-step*(tr-1), 'Label',sprintf('Trial%d',tr));
            yl.LabelHorizontalAlignment = 'center';
            yl.LabelVerticalAlignment = 'bottom';
        end
        xlabel('Time (ms)');
        ylabel('Amplitude (mV)');
        title(sprintf('%s - %s', subject, EEG_epoch2.chanlocs(:).labels))
    end
end