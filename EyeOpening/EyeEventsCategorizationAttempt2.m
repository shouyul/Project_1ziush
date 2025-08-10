clear all;
config_PIONEER_Tumbler;

eyeDetectionChannel = {'R3Z'};
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
    subject_ind = 4;
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
    EEG_epoch1 = pop_epoch(EEG_filt1, {'OpenEyesSoundOn'}, [-7.5,27.5]);
    
    EEG_filt2 = custom_filter(EEG_sel, 0.05, 5);
    EEG_epoch2 = pop_epoch(EEG_filt2, {'OpenEyesSoundOn'}, [-7.5,27.5]);
    
    % General parameters
    opts_peaks.timeEO = 0;
    opts_peaks.timeEC = 20860;
    opts_peaks.accept_wdw = [-1,5]; % in seconds
    opts_peaks.percMax4Prom = 0;
    opts_peaks.wdw_diff1 = 0.2; % in seconds, for diff calculation
    opts_peaks.wdw_diff2 = 0.4; % in seconds, for diff calculation
    opts_peaks.diffQuant = 75; % for createPeaksTable2
    opts_peaks.artifact_detection_th = 500; % in microV
    
    % Train set
    n_trials_train = 10;
    train_trials = randi([1 EEG_epoch1.trials],1,n_trials_train);
    while length(train_trials) ~= length(unique(train_trials))
        train_trials = randi([1 EEG_epoch1.trials],1,n_trials_train);
    end
    
    opts_peaks.quant = 90;
    opts_peaks.userLabel = true;
    opts_peaks.inspectedTrials = train_trials;
    Peaks_train = createPeaksTable(EEG_epoch1, EEG_epoch2, train_trials, opts_peaks);
    Peaks = createPeaksTable2(EEG_epoch1, EEG_epoch2, opts_peaks);
    
    % Test set
    %test_trials = 1:EEG_epoch1.trials;
    test_trials = setdiff(1:EEG_epoch1.trials,train_trials);
    opts_peaks.quant = 95;
    opts_peaks.userLabel = false;
    Peaks_test = createPeaksTable(EEG_epoch1, EEG_epoch2, test_trials, opts_peaks);
    
    % Classify
    cats_binary = cell(size(Peaks_train.EyeOpening));
    cats_binary(Peaks_train.Rejected) = {'NoEvent'};
    cats_binary(~Peaks_train.Rejected) = {'Event'};
    
    cats_full = cell(size(Peaks_train.EyeOpening));
    cats_full(Peaks_train.Rejected) = {'NoEvent'};
    cats_full(Peaks_train.EyeOpening) = {'EyeOpening'};
    cats_full(Peaks_train.EyeClosing) = {'EyeClosing'};
    cats_full(Peaks_train.Blink) = {'Blink'};
      
    Peaks_train = addvars(Peaks_train, categorical(cats_binary), categorical(cats_full),...
        'NewVariableNames', {'Category', 'SubCategory'});    
    Peak_train_events = Peaks_train(~Peaks_train.Rejected,:);    
    
    formula1 = 'Category~Latency+Height1+Width1+Prominence1+Diff1+Diff2';
    formula2 = 'SubCategory~Latency+Height1+Width1+Prominence1+Diff1+Diff2';
    cost_bin = ones(2) - diag([1 1]);
    cost_full = ones(4) - diag([1 1 1 1]);
    cost_full(:,1) = cost_full(:,1)*10; % Predicting rejection for other classes
    cost_full(1,:) = cost_full(1,:)/10; % Predicting an other class for a true rejection
    
    simpleTree1 = fitctree(Peaks_train, formula1, 'Cost', cost_bin);
    view(simpleTree1,'Mode', 'graph')
    simpleTree2 = fitctree(Peaks_train, formula2, 'Cost', cost_full);
    view(simpleTree2,'Mode', 'graph')
    
    optimizedTree1 = fitctree(Peaks_train, formula1, 'Cost', cost_bin,...
        'OptimizeHyperparameters','auto');
    classError1 = loss(optimizedTree1, Peaks_train, 'Category');
    view(optimizedTree1,'Mode', 'graph')
    optimizedTree2 = fitctree(Peaks_train, formula2, 'Cost', cost_full,...
        'OptimizeHyperparameters','auto');
    classError2 = loss(optimizedTree2, Peaks_train, 'SubCategory');
    view(optimizedTree2,'Mode', 'graph')
    
    % Predict
    labels = predict(optimizedTree2, Peaks_test);
    Peaks_test.Rejected = labels == 'NoEvent';
    Peaks_test.EyeOpening = labels == 'EyeOpening';
    Peaks_test.EyeClosing = labels == 'EyeClosing';
    Peaks_test.Blink = labels == 'Blink';
    
    % Plot
    opts_plot.subject = subject;
    opts_plot.step = 0.25; % in mV
    opts_plot.timeEO = 0;
    opts_plot.timeEC = 20860;
    %plotEyeEventsOnTrials(EEG_epoch1, EEG_epoch2, Peaks_test, opts_plot)
    plotEyeEventsOnTrials(EEG_epoch1, EEG_epoch2, Peaks, opts_plot)
    
    continue
    
    % General parameters
    n_trials_train = 10;
    artifact_thresh = 500; % in microV
    quant = 90;
    percMax4Prom = 0;
    wdw_diff1 = 0.2; % in seconds, for diff calculation
        latsLeft1 = -round(wdw_diff1*EEG.srate):-1;
    latsRight1 = 1:round(wdw_diff1*EEG.srate);
    wdw_diff2 = 0.4; % in seconds, for diff calculation
            latsLeft2 = -round(wdw_diff1*EEG.srate):-1;
    latsRight2 = 1:round(wdw_diff1*EEG.srate);
    diffminEO = 5;
    diffmaxEC = -10;
    
    wdw_lat_bef = 1; % in seconds, for latency criteria
    wdw_lat_aft = 5; % in seconds, for latency criteria
    refLatEO = find(EEG_epoch2.times == 0);
    refLatEC = find(EEG_epoch2.times == 20860);
    
    varNames = {'Trial','Quantile','Latency','Max1','Ave1',...
        'Height1','Width1','Prominence1','Max2','Ave2','Height2','Diff1','Diff2',...
        'EyeOpening','EyeClosing','Blink','Rejected'};
    varTypes = {'int32','double','int32','double','double',...
        'double','double','double','double','double','double','double','double',...
        'logical','logical','logical','logical'};
    
    Peaks_train = table('Size',[0,numel(varNames)], 'VariableTypes',varTypes, 'VariableNames', varNames);
    i_train = 0;
    
    train_trials = randi([1 EEG_epoch1.trials],1,n_trials_train);
    while length(train_trials) ~= length(unique(train_trials))
        train_trials = randi([1 EEG_epoch1.trials],1,n_trials_train);
    end
    
    % Train
    for tr = train_trials
        EEGave1 = mean(EEG_epoch1.data(:,:,tr), 'all');
        EEGmax1 = max(EEG_epoch1.data(:,:,tr), [], 'all');
        EEGquant = quantile(EEG_epoch1.data(:,:,tr), quant/100, 'all');
        [heights, locs, widths, proms] = findpeaks(squeeze(EEG_epoch1.data(:,:,tr)),...
            'MinPeakHeight', EEGquant, 'MinPeakProminence', EEGmax1*percMax4Prom);
        
        EEGave2 = mean(EEG_epoch2.data(:,:,tr), 'all');
        EEGmax2 = max(EEG_epoch2.data(:,:,tr), [], 'all');
        
        Npeaks = length(locs);
        for p = 1:Npeaks
            i_train = i_train+1;
            Peaks_train.Trial(i_train) = tr;
            Peaks_train.Quantile(i_train) = quant;
            Peaks_train.Latency(i_train) = locs(p);
            Peaks_train.Max1(i_train) = EEGmax1;
            Peaks_train.Ave1(i_train) = EEGave1;
            Peaks_train.Height1(i_train) = heights(p);
            Peaks_train.Width1(i_train) = widths(p);
            Peaks_train.Prominence1(i_train) = proms(p);
            Peaks_train.Max2(i_train) = EEGmax2;
            Peaks_train.Ave2(i_train) = EEGave2;
            
            % Default values
            Peaks_train.EyeOpening(i_train) = false;
            Peaks_train.EyeClosing(i_train) = false;
            Peaks_train.Blink(i_train) = false;
            Peaks_train.Rejected(i_train) = false;
            if locs(p) > wdw_diff2*EEG.srate && locs(p) < EEG_epoch1.pnts - wdw_diff2*EEG.srate
                Peaks_train.Height2(i_train) = mean(EEG_epoch2.data(:,locs(p),tr),'all');
                
                %                 if Peaks_train.Height2(i_train) >= artifact_thresh
                %                     Peaks_train.Diff2(i_train) = NaN;
                %                     Peaks_train.Rejected(i_train) = true;
                %                     continue
                %                 end
                
                meanLeft1 = mean(EEG_epoch2.data(:,locs(p)+latsLeft1,tr),'all');
                meanRight1 = mean(EEG_epoch2.data(:,locs(p)+latsRight1,tr),'all');
                diffLR1 = meanLeft1-meanRight1;
                Peaks_train.Diff1(i_train) = diffLR1;
                
                                meanLeft2 = mean(EEG_epoch2.data(:,locs(p)+latsLeft2,tr),'all');
                meanRight2 = mean(EEG_epoch2.data(:,locs(p)+latsRight2,tr),'all');
                diffLR2 = meanLeft2-meanRight2;
                Peaks_train.Diff2(i_train) = diffLR2;
                
                %                 if diffLR > diffminEO
                %                     Peaks_train.EyeOpening(i_train) = true;
                %                 elseif diffLR < diffmaxEC
                %                     Peaks_train.EyeClosing(i_train) = true;
                %                 else
                %                     Peaks_train.Blink(i_train) = true;
                %                 end
            else
                Peaks_train.Height2(i_train) = NaN;
                Peaks_train.Diff1(i_train) = NaN;
                Peaks_train.Diff2(i_train) = NaN;
                Peaks_train.Rejected(i_train) = true;
            end
        end
        
        figure;hold on;
        plot(EEG_epoch1.times,mean(EEG_epoch1.data(:,:,tr),1),':k');
        plot(EEG_epoch2.times,mean(EEG_epoch2.data(:,:,tr),1));
        xline(0, 'Label', 'Eye Opening event');
        xline(20860, 'Label', 'Eye Closing event');
        yline(0,'--');
        yline(EEGquant,':');
        xlabel('Time (ms)');
        ylabel('Amplitude (microV)');
        title(sprintf('Trial %d', tr));
        
        peaks_tr = find(Peaks_train.Trial == tr)';
        scatter(EEG_epoch1.times(Peaks_train.Latency(peaks_tr)), Peaks_train.Height1(peaks_tr),'k','*');
        for p = peaks_tr
            if ~isnan(Peaks_train.Height2(p))
                pt = scatter(EEG_epoch2.times(Peaks_train.Latency(p)), Peaks_train.Height2(p),'m','*');
                accepted = false;
                while ~accepted
                    dec = input('Which category? [0=Exclude; 1=EyeOpening; 2=EyeClosing; 3=Blink] ');
                    if ~isempty(dec)
                        switch dec
                            case 0
                                accepted = true;
                                Peaks_train.Rejected(p) = true;
                                delete(pt);
                                scatter(EEG_epoch2.times(Peaks_train.Latency(p)), Peaks_train.Height2(p),'r','+');
                            case 1
                                accepted = true;
                                Peaks_train.EyeOpening(p) = true;
                                delete(pt);
                                scatter(EEG_epoch2.times(Peaks_train.Latency(p)), Peaks_train.Height2(p),'g','o');
                            case 2
                                accepted = true;
                                Peaks_train.EyeClosing(p) = true;
                                delete(pt);
                                scatter(EEG_epoch2.times(Peaks_train.Latency(p)), Peaks_train.Height2(p),'g','x');
                            case 3
                                accepted = true;
                                Peaks_train.Blink(p) = true;
                                delete(pt);
                                scatter(EEG_epoch2.times(Peaks_train.Latency(p)), Peaks_train.Height2(p),'g','v');
                            otherwise
                                fprintf('Decision not understood.\n')
                        end
                    end
                end
            end
        end
        %close gcf
        
        %         % Look at the sequencing
        %         tr_peaks = Peaks_train.Trial == tr;
        %         eo_peaks = find(Peaks_train.EyeOpening & tr_peaks);
        %         ec_peaks = find(Peaks_train.EyeClosing & tr_peaks);
        %
        %         Peaks_train.Rejected(eo_peaks) = true;
        %         Peaks_train.Rejected(ec_peaks) = true;
        %         [diff_eo,sortorder_eo] = sort(Peaks_train.Diff2(eo_peaks),'descend');
        %         [diff_ec,sortorder_ec] = sort(Peaks_train.Diff2(ec_peaks),'ascend');
        %
        %         lats = 1:EEG_epoch1.pnts;
        %         eo_period = false(1,length(lats));
        %         c_max = max([length(eo_peaks),length(ec_peaks)]);
        %         for c = 1:c_max
        %             if c == 1
        %                 if length(eo_peaks) >= 1 && length(ec_peaks) >= 1
        %                     Peaks_train.Rejected(eo_peaks(sortorder_eo(c))) = false;
        %                     Peaks_train.Rejected(ec_peaks(sortorder_ec(c))) = false;
        %                     if Peaks_train.Latency(eo_peaks(sortorder_eo(c))) < Peaks_train.Latency(ec_peaks(sortorder_ec(c)))
        %                         eo_period(lats >= Peaks_train.Latency(eo_peaks(sortorder_eo(c))) &...
        %                             lats < Peaks_train.Latency(ec_peaks(sortorder_ec(c)))) = true;
        %                     else
        %                         eo_period(lats >= Peaks_train.Latency(eo_peaks(sortorder_eo(c))) |...
        %                             lats < Peaks_train.Latency(ec_peaks(sortorder_ec(c)))) = true;
        %                     end
        %                 elseif length(eo_peaks) >= 1
        %                     Peaks_train.Rejected(eo_peaks(sortorder_eo(c))) = false;
        %                     eo_period(lats >= Peaks_train.Latency(eo_peaks(sortorder_eo(c)))) = true;
        %                 elseif length(ec_peaks) >= 1
        %                     Peaks_train.Rejected(ec_peaks(sortorder_ec(c))) = false;
        %                     eo_period(lats < Peaks_train.Latency(ec_peaks(sortorder_ec(c)))) = true;
        %                 else
        %                     eo_period(lats) = true;
        %                 end
        %             else
        %                 break
        %                 % Not completed yet
        %                 if c <= length(eo_peaks) && ec_period(Peaks_train.Latency(eo_peaks(sortorder_eo(c))))
        %                     Peaks_train.Rejected(eo_peaks(sortorder_eo(c))) = false;
        %
        %                 end
        %             end
        %         end
        
        %         [~,best_eo] = max(Peaks.Diff2(eo_peaks));
        %         Peaks.Rejected(eo_peaks(best_eo)) = false;
        %         [~,best_ec] = min(Peaks.Diff2(ec_peaks));
        %         Peaks.Rejected(ec_peaks(best_ec)) = false;
        
        %         bl_peaks = find(Peaks_train.Blink & tr_peaks);
        %         Peaks_train.Rejected(bl_peaks(~eo_period(Peaks_train.Latency(bl_peaks)))) = true;
        %
        %         if Peaks.Latency(eo_peaks(best_eo)) < Peaks.Latency(ec_peaks(best_ec))
        %             Peaks.Rejected(bl_peaks(Peaks.Latency(bl_peaks) < Peaks.Latency(eo_peaks(best_eo)) |...
        %                 Peaks.Latency(bl_peaks) > Peaks.Latency(ec_peaks(best_ec)))) = true;
        %         else
        %             Peaks.Rejected(bl_peaks(Peaks.Latency(bl_peaks) < Peaks.Latency(eo_peaks(best_eo)) &...
        %                 Peaks.Latency(bl_peaks) > Peaks.Latency(ec_peaks(best_ec)))) = true;
        %         end
        
        % Plotting
        
        
        %             scatter(EEG_epoch2.times(Peaks_train.Latency(Peaks_train.EyeOpening & ~Peaks_train.Rejected & tr_peaks)),...
        %                 Peaks_train.Height2(Peaks_train.EyeOpening & ~Peaks_train.Rejected & tr_peaks),'g','o')
        %             scatter(EEG_epoch2.times(Peaks_train.Latency(Peaks_train.EyeOpening & Peaks_train.Rejected & tr_peaks)),...
        %                 Peaks_train.Height2(Peaks_train.EyeOpening & Peaks_train.Rejected & tr_peaks),'r','o')
        %             scatter(EEG_epoch2.times(Peaks_train.Latency(Peaks_train.EyeClosing & ~Peaks_train.Rejected & tr_peaks)),...
        %                 Peaks_train.Height2(Peaks_train.EyeClosing & ~Peaks_train.Rejected & tr_peaks),'g','x')
        %             scatter(EEG_epoch2.times(Peaks_train.Latency(Peaks_train.EyeClosing & Peaks_train.Rejected & tr_peaks)),...
        %                 Peaks_train.Height2(Peaks_train.EyeClosing & Peaks_train.Rejected & tr_peaks),'r','x')
        %             scatter(EEG_epoch2.times(Peaks_train.Latency(Peaks_train.Blink & ~Peaks_train.Rejected & tr_peaks)),...
        %                 Peaks_train.Height2(Peaks_train.Blink & ~Peaks_train.Rejected & tr_peaks),'g','*')
        %             scatter(EEG_epoch2.times(Peaks_train.Latency(Peaks_train.Blink & Peaks_train.Rejected & tr_peaks)),...
        %                 Peaks_train.Height2(Peaks_train.Blink & Peaks_train.Rejected & tr_peaks),'r','*')
        %             xline(0, 'Label', 'Eye Opening event');
        %             xline(20860, 'Label', 'Eye Closing event');
        %             yline(0,'--');
        %             yline(EEGquant,':');
        %             xlabel('Time (ms)');
        %             ylabel('Amplitude (microV)');
        %             title(sprintf('Trial %d', tr));
        
    end
    
    cats_full = double(Peaks_train.EyeOpening);
    cats_full(Peaks_train.EyeClosing) = 2;
    cats_full(Peaks_train.Blink) = 3;
    cats_full(Peaks_train.Rejected) = 0;    
    Peaks_train = addvars(Peaks_train, categorical(cats_full), 'NewVariableNames', {'Category'});
    %categories_SB = categories_full > 0;
    formula1 = 'Category~Latency+Height1+Width1+Prominence1+Height2+Diff1+Diff2';
    cost_full = ones(4) - diag([1 1 1 1]);
    cost_full(:,1) = cost_full(:,1)*10; % Predicting rejection for other classes
    cost_full(1,:) = cost_full(1,:)/10; % Predicting an other class for a true rejection
    simpleTree = fitctree(Peaks_train, formula1, 'Cost', cost_full);
    view(simpleTree,'Mode', 'graph')
    optimizedTree = fitctree(Peaks_train, formula1, 'Cost', cost_full,...
        'OptimizeHyperparameters','auto');
    view(optimizedTree,'Mode', 'graph')
    
%     tree = fitctree(Peaks_train,'Category~Latency+Height1+Width1+Prominence1+Height2+Diff2',...
%         'OptimizeHyperparameters','auto');
%     view(tree.Trained{1}, 'Mode', 'graph')
%     view(tree.Trained{2}, 'Mode', 'graph')
%     view(tree.Trained{3}, 'Mode', 'graph')
%     view(tree,'Mode', 'graph')
%     
%     tree1 = prune(tree);
    

Peaks_test = table('Size',[0,numel(varNames)], 'VariableTypes',varTypes, 'VariableNames', varNames);
    i_test = 0;
    
    test_trials = setdiff(1:EEG_epoch1.trials,train_trials);
    
    % Test
    for tr = test_trials
        EEGave1 = mean(EEG_epoch1.data(:,:,tr), 'all');
        EEGmax1 = max(EEG_epoch1.data(:,:,tr), [], 'all');
        EEGquant = quantile(EEG_epoch1.data(:,:,tr), quant/100, 'all');
        [heights, locs, widths, proms] = findpeaks(squeeze(EEG_epoch1.data(:,:,tr)),...
            'MinPeakHeight', EEGquant, 'MinPeakProminence', EEGmax1*percMax4Prom);
        
        EEGave2 = mean(EEG_epoch2.data(:,:,tr), 'all');
        EEGmax2 = max(EEG_epoch2.data(:,:,tr), [], 'all');
        
        Npeaks = length(locs);
        for p = 1:Npeaks
            i_test = i_test+1;
            Peaks_test.Trial(i_test) = tr;
            Peaks_test.Quantile(i_test) = quant;
            Peaks_test.Latency(i_test) = locs(p);
            Peaks_test.Max1(i_test) = EEGmax1;
            Peaks_test.Ave1(i_test) = EEGave1;
            Peaks_test.Height1(i_test) = heights(p);
            Peaks_test.Width1(i_test) = widths(p);
            Peaks_test.Prominence1(i_test) = proms(p);
            Peaks_test.Max2(i_test) = EEGmax2;
            Peaks_test.Ave2(i_test) = EEGave2;
            
            % Default values
            Peaks_test.EyeOpening(i_test) = false;
            Peaks_test.EyeClosing(i_test) = false;
            Peaks_test.Blink(i_test) = false;
            Peaks_test.Rejected(i_test) = false;
            if locs(p) > wdw_diff2*EEG.srate && locs(p) < EEG_epoch1.pnts - wdw_diff2*EEG.srate
                Peaks_test.Height2(i_test) = mean(EEG_epoch2.data(:,locs(p),tr),'all');
                
                %                 if Peaks_train.Height2(i_train) >= artifact_thresh
                %                     Peaks_train.Diff2(i_train) = NaN;
                %                     Peaks_train.Rejected(i_train) = true;
                %                     continue
                %                 end
                
                
                                meanLeft1 = mean(EEG_epoch2.data(:,locs(p)+latsLeft1,tr),'all');
                meanRight1 = mean(EEG_epoch2.data(:,locs(p)+latsRight1,tr),'all');
                diffLR1 = meanLeft1-meanRight1;
                Peaks_test.Diff1(i_test) = diffLR1;
                
                                meanLeft2 = mean(EEG_epoch2.data(:,locs(p)+latsLeft2,tr),'all');
                meanRight2 = mean(EEG_epoch2.data(:,locs(p)+latsRight2,tr),'all');
                diffLR2 = meanLeft2-meanRight2;
                Peaks_test.Diff2(i_test) = diffLR2;
                
                %                 if diffLR > diffminEO
                %                     Peaks_train.EyeOpening(i_train) = true;
                %                 elseif diffLR < diffmaxEC
                %                     Peaks_train.EyeClosing(i_train) = true;
                %                 else
                %                     Peaks_train.Blink(i_train) = true;
                %                 end
            else
                Peaks_test.Height2(i_test) = NaN;
                Peaks_test.Diff1(i_test) = NaN;
                Peaks_test.Diff2(i_test) = NaN;
                Peaks_test.Rejected(i_test) = true;
            end
        end
    end
    
    labels = predict(optimizedTree, Peaks_test);
    Peaks_test.Rejected = labels == '0';
    Peaks_test.EyeOpening = labels == '1';
    Peaks_test.EyeClosing = labels == '2';
    Peaks_test.Blink = labels == '3';


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
        edges = 0:5:max(Peaks_train.Height1)+5;
        histogram(Peaks_train.Height1(Peaks_train.EyeOpening & ~Peaks_train.Rejected),edges);
        histogram(Peaks_train.Height1(Peaks_train.EyeClosing & ~Peaks_train.Rejected),edges);
        histogram(Peaks_train.Height1(Peaks_train.Blink & ~Peaks_train.Rejected),edges);
        %histogram(Peaks.Height1(Peaks.Rejected),edges);
        title('Height1')
        subplot(3,2,3)
        hold on;
        edges = 0:10:max(Peaks_train.Width1)+10;
        histogram(Peaks_train.Width1(Peaks_train.EyeOpening & ~Peaks_train.Rejected),edges);
        histogram(Peaks_train.Width1(Peaks_train.EyeClosing & ~Peaks_train.Rejected),edges);
        histogram(Peaks_train.Width1(Peaks_train.Blink & ~Peaks_train.Rejected),edges);
        %histogram(Peaks.Width1(Peaks.Rejected),edges);
        title('Width1')
        subplot(3,2,5)
        hold on;
        edges = 0:5:max(Peaks_train.Prominence1)+5;
        histogram(Peaks_train.Prominence1(Peaks_train.EyeOpening & ~Peaks_train.Rejected),edges);
        histogram(Peaks_train.Prominence1(Peaks_train.EyeClosing & ~Peaks_train.Rejected),edges);
        histogram(Peaks_train.Prominence1(Peaks_train.Blink & ~Peaks_train.Rejected),edges);
        %histogram(Peaks.Prominence1(Peaks.Rejected),edges);
        title('Prominence1')
        subplot(3,2,2)
        hold on;
        edges = 0:10:max(Peaks_train.Height2)+10;
        histogram(Peaks_train.Height2(Peaks_train.EyeOpening & ~Peaks_train.Rejected),edges);
        histogram(Peaks_train.Height2(Peaks_train.EyeClosing & ~Peaks_train.Rejected),edges);
        histogram(Peaks_train.Height2(Peaks_train.Blink & ~Peaks_train.Rejected),edges);
        %histogram(Peaks.Height2(Peaks.Rejected),edges);
        title('Height2')
        subplot(3,2,4)
        hold on;
        edges = (min(Peaks_train.Diff2)-10):10:max(Peaks_train.Diff2)+10;
        h1 = histogram(Peaks_train.Diff2(Peaks_train.EyeOpening & ~Peaks_train.Rejected),edges);
        h2 = histogram(Peaks_train.Diff2(Peaks_train.EyeClosing & ~Peaks_train.Rejected),edges);
        h3 = histogram(Peaks_train.Diff2(Peaks_train.Blink & ~Peaks_train.Rejected),edges);
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
        
        for tr = test_trials
            decMeanData1(:,tr) = meanData1(:,tr) - step*(tr-1);
            decMeanData2(:,tr) = meanData2(:,tr) - step*(tr-1);
        end
        
        peaks4plot = nan(size(Peaks_test,1),1);
        for p = 1:size(Peaks_test,1)
            tr = Peaks_test.Trial(p);
            peaks4plot(p) = decMeanData2(Peaks_test.Latency(p),tr);
        end
        
        figure;hold on;
        plot(EEG_epoch1.times,decMeanData1,':k');
        plot(EEG_epoch2.times,decMeanData2);
        scatter(EEG_epoch2.times(Peaks_test.Latency(Peaks_test.EyeOpening)),...
            peaks4plot(Peaks_test.EyeOpening),'g','o')
        scatter(EEG_epoch2.times(Peaks_test.Latency(Peaks_test.EyeClosing)),...
            peaks4plot(Peaks_test.EyeClosing),'g','x')
        scatter(EEG_epoch2.times(Peaks_test.Latency(Peaks_test.Blink)),...
            peaks4plot(Peaks_test.Blink),'g','v')
        scatter(EEG_epoch2.times(Peaks_test.Latency(Peaks_test.Rejected)),...
            peaks4plot(Peaks_test.Rejected),'r','*')
        
        xline(0, 'Label', 'Eye Opening event');
        xline(20850, 'Label', 'Eye Closing event');
        for tr = test_trials
            yl = yline(-step*(tr-1), 'Label',sprintf('Trial%d',tr));
            yl.LabelHorizontalAlignment = 'center';
            yl.LabelVerticalAlignment = 'bottom';
        end
        xlabel('Time (ms)');
        ylabel('Amplitude (mV)');
        title(sprintf('%s - %s', subject, EEG_epoch2.chanlocs(:).labels))
    end
end