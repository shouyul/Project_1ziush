clear all;
config_PIONEER_Tumbler;

for subject_ind = subject_inds
    if ~exist('EEG', 'var')
        switch user
            case 'sl'
                launchEEGLAB
            case 'JB'
                eeglab
        end
    end
    
    subject = study_config.subjects(subject_ind).id;
    disp(['Subject ' subject]);
    study_config.current_subject = subject_ind;
    
    %% Folders and Files names:
    N = makeFolderFileNames(study_config, subject);
    input_filepath = N.searchFolder_3arch_rej_ICcats;
    
    if strcmp(study_config.epochs.event, 'TumblerVisible')
        EEG = pop_loadset('filename', N.epochedFile, 'filepath', input_filepath);
        chans2select = {'Z17Z','L18Z','R18Z','Z19Z'};
        % Select latencies
        EEG = pop_select(EEG, 'time', study_config.epochs.limits_wdw);
        
        TumbVisChan = strcmp({EEG.chanlocs.labels}, 'TumblerVisibility');
        figure;
        subplot(1+numel(chans2select),1,1)
        plot(EEG.times,squeeze(EEG.data(TumbVisChan,:,:)));
        ylim([0,100])
        xline(0,'--k');
        title('Tumbler Visibility')
        
        % Select channels
        EEG = pop_select(EEG, 'channel', chans2select);
        
        % options for baseline comparison
        opts_stats.model = 'classic';
        opts_stats.pairing = 'on';
        opts_stats.style = 'time';
        opts_stats.N_reps = 1000;
        opts_stats.reusePerms = false;
        opts_stats.removeSmallestClusters = false;
        opts_stats.fields = {'Visible','Visible_base'};
        %         stat_params.apply_clusterCorr = false;
        %         stat_params.correctionType = 'bonf';
        %         stat_params.apply_pairwiseCorr = true;
        
        %ERPs_stats.Subjects = subjects;
        ERPs_stats.Times = EEG.times;
        ERPs_stats.BaselineModel = 'None';
        
        base_times = EEG.times<0;
        
        for ch = 1:EEG.nbchan
            fprintf('Computing statistics vs baseline\n')
            data4stats = squeeze(EEG.data(ch,:,:));
            
            ERPs_base_surrog = nan(size(data4stats));
            % Create baseline surrogate trials
            base_distrib = data4stats(base_times,:);
            base_distrib = base_distrib(:);
            % Simply draw some samples from the distribution
            draw = ceil(length(base_distrib)*...
                rand(size(ERPs_base_surrog)));
            ERPs_base_surrog = base_distrib(draw);
            
            ERPs_stats.Visible = data4stats;
            ERPs_stats.Visible_base = ERPs_base_surrog;
            
            masks = struct();
            masks.n_masks = 1;
            masks.nameMask1 = 'Clustered T-stats Masks';
            
            %%%%% Clustered statistics %%%%%
            [clustered_stats_table, statistical_clusters, stats_surrog, ~, permutations] =...
                NP_statTest(ERPs_stats, opts_stats);
            
            signif_clusters = clustered_stats_table.ClusterID(clustered_stats_table.ClusterPval <= 0.05);
            clustered_mask = false(size(statistical_clusters));
            for signCl = signif_clusters'
                clustered_mask = clustered_mask | statistical_clusters == signCl;
            end
            masks.mask1 = clustered_mask;
            
            operations = struct();
            operations.singleTrialNorm = false;
            operations.preStimBaseline = false;
            operations.baseline = [];
            opts_plot = struct();
            opts_plot.style = 'GrandAverage';
            opts_plot.Smooth = false;
            opts_plot.nSampsSmooth = 25;
            opts_plot.ylims = [-7.5,7.5];
            if ch == EEG.nbchan
            opts_plot.xlabel = true;
            else
                opts_plot.xlabel = false;
            end
            opts_plot.ylabel = true;
            
            subplot(1+EEG.nbchan,1,1+ch)
            plotERP(EEG.data(ch,:,:), EEG.times, masks.mask1, operations, opts_plot);
            title(sprintf('ERP on %s', EEG.chanlocs(ch).labels))
            suptitle(subject);
        end
        
        continue
        
        mocapChans = contains({EEG.chanlocs.type}, 'MOCAP');
        nb_chans = EEG.nbchan - sum(mocapChans);
        %         if any(mocapChans)
        %             nb_chans = EEG.nbchan - sum(mocapChans);
        %             EEG = pop_select(EEG, 'nochannel', mocapChans);
        %         else
        %             nb_chans = EEG.nbchan;
        %         end
        
        nb_trials = EEG.trials;
        freqs = study_config.psd.FoI;
        nb_freqs = length(freqs);
        %SR = EEG.srate;
        
        TrialsInfo = EEG.etc.epochInfo.Trials;
        TrialsInfo = TrialsInfo(:,1:8);
        
        [PSD_EC, TumbVis_EC, nb_chunks_EC] = computePSD_TumbVis(EEG, 'EC', study_config);
        [PSD_EO, TumbVis_EO, nb_chunks_EO] = computePSD_TumbVis(EEG, 'EO', study_config);
        
        %         if study_config.psd.chunks == 0
        %             nb_chunks_EC = 1;
        %             nb_chunks_EO = 1;
        %         else
        %             nb_chunks_EC = abs(study_config.epochs.limits_wdw(1))/study_config.psd.chunks;
        %             nb_chunks_EO = study_config.epochs.limits_wdw(2)/study_config.psd.chunks;
        %         end
        
        %         PSD_EC = zeros(nb_chans, nb_freqs, nb_trials, nb_chunks_EC);
        %         PSD_EO = zeros(nb_chans, nb_freqs, nb_trials, nb_chunks_EO);
        %         if any(mocapChans)
        %             TumbVisChan = strcmp({EEG.chanlocs.labels}, 'TumblerVisibility');
        %             TumbVis_EC = zeros(nb_trials, nb_chunks_EC);
        %             TumbVis_EO = zeros(nb_trials, nb_chunks_EO);
        %         end
        
        %         switch study_config.psd.method
        %             case 'pwelch'
        %                 % Compute PSD with pwelch method
        %                 overlap = SR/25;
        %                 if nb_chunks_EC == 1
        %                     for tr = 1:nb_trials
        %                         PSD_EC(:,:,tr,1) = pwelch(squeeze(EEG_EC.data(~mocapChans,:,tr))', SR, overlap, freqs, SR)';
        %
        %                         if any(mocapChans)
        %                             TumbVis_EC(tr,1) = mean(squeeze(EEG_EC.data(TumbVisChan,:,tr)));
        %                         end
        %                     end
        %                 else
        %                     for ch = 1:nb_chunks_EC
        %                         EEG_EC_sel = pop_select(EEG, 'time', study_config.epochs.limits_wdw(1)+...
        %                             [ch-1,ch].*study_config.psd.chunks);
        %                         for tr = 1:nb_trials
        %                             PSD_EC(:,:,tr,ch) = pwelch(squeeze(EEG_EC_sel.data(~mocapChans,:,tr))', SR, overlap, freqs, SR)';
        %                             if any(mocapChans)
        %                                 TumbVis_EC(tr,ch) = mean(squeeze(EEG_EC_sel.data(TumbVisChan,:,tr)));
        %                             end
        %                         end
        %                     end
        %                 end
        %
        %                 if nb_chunks_EO == 1
        %                     for tr = 1:nb_trials
        %                         PSD_EO(:,:,tr,1) = pwelch(squeeze(EEG_EO.data(~mocapChans,:,tr))', SR, overlap, freqs, SR)';
        %                         if any(mocapChans)
        %                             TumbVis_EO(tr,1) = mean(squeeze(EEG_EO.data(TumbVisChan,:,tr)));
        %                         end
        %                     end
        %                 else
        %                     for ch = 1:nb_chunks_EO
        %                         EEG_EO_sel = pop_select(EEG, 'time', [ch-1,ch].*study_config.psd.chunks);
        %                         for tr = 1:nb_trials
        %                             PSD_EO(:,:,tr,ch) = pwelch(squeeze(EEG_EO_sel.data(~mocapChans,:,tr))', SR, overlap, freqs, SR)';
        %                             if any(mocapChans)
        %                                 TumbVis_EO(tr,ch) = mean(squeeze(EEG_EO_sel.data(TumbVisChan,:,tr)));
        %                             end
        %                         end
        %                     end
        %                 end
        %
        %             otherwise
        %                 error('Not coded yet')
        %         end
        
        TrialsInfo_EC = adaptTrialsInfo(TrialsInfo, TumbVis_EC, nb_chunks_EC);
        PSD_EC = reshape(permute(PSD_EC,[1,2,4,3]),[nb_chans, nb_freqs, nb_trials*nb_chunks_EC]);
        TrialsInfo_EO = adaptTrialsInfo(TrialsInfo, TumbVis_EO, nb_chunks_EO);
        PSD_EO = reshape(permute(PSD_EO,[1,2,4,3]),[nb_chans, nb_freqs, nb_trials*nb_chunks_EO]);
        
        % Concatenate over subjects:
        TrialsInfo_EC_all = cat(1, TrialsInfo_EC_all, TrialsInfo_EC);
        PSD_EC_all = cat(3,PSD_EC_all,PSD_EC);
        TrialsInfo_EO_all = cat(1, TrialsInfo_EO_all, TrialsInfo_EO);
        PSD_EO_all = cat(3,PSD_EO_all,PSD_EO);
        
        %% Normalize data according to options
        [PSD_EC_norm] = normalizePSD(PSD_EC, TrialsInfo_EC, study_config.psd);
        [PSD_EO_norm] = normalizePSD(PSD_EO, TrialsInfo_EO, study_config.psd);
        
        %% Save data
        opts_plot = study_config.psd;
        opts_plot.SR = EEG.srate;
        folderName = makePSDFolderName(opts_plot);
        if ~exist(fullfile(input_filepath,folderName),'dir')
            mkdir(fullfile(input_filepath,folderName));
        end
        opts_plot.event = study_config.epochs.event;
        opts_plot.phase = 'EC';
        
        fileName = makePSDFileName('labels', subject, opts_plot);
        save(fullfile(input_filepath, fileName), sprintf('TrialsInfo_%s',opts_plot.phase));
        
        fileName = makePSDFileName('psd', subject, opts_plot);
        save(fullfile(input_filepath, folderName, fileName), sprintf('PSD_%s_norm',opts_plot.phase), '-v7.3');
        
        opts_plot.phase = 'EO';
        fileName = makePSDFileName('labels', subject, opts_plot);
        save(fullfile(input_filepath, fileName), sprintf('TrialsInfo_%s',opts_plot.phase));
        
        fileName = makePSDFileName('psd', subject, opts_plot);
        save(fullfile(input_filepath, folderName, fileName), sprintf('PSD_%s_norm',opts_plot.phase), '-v7.3');
    else
        error('Event not suitable for this analysis')
    end
end

%% Multi-subjects data
% if strcmp(study_config.epochs.event, 'EyesOpening')
%     %% Normalize data according to options
%     [PSD_EC_all_norm] = normalizePSD(PSD_EC_all, TrialsInfo_EC_all, study_config.psd);
%     [PSD_EO_all_norm] = normalizePSD(PSD_EO_all, TrialsInfo_EO_all, study_config.psd);
%
%     %% Save data
%     if ~exist(fullfile(N.searchFolder_4arch_rej_ICcats,folderName),'dir')
%         mkdir(fullfile(N.searchFolder_4arch_rej_ICcats,folderName));
%     end
%     options.phase = 'EC';
%     fileName = makePSDFileName('labels', 'AllSubjs', options);
%     save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), sprintf('TrialsInfo_%s_all',options.phase));
%
%     fileName = makePSDFileName('psd', 'AllSubjs', options);
%     save(fullfile(N.searchFolder_4arch_rej_ICcats, folderName, fileName), sprintf('PSD_%s_all_norm', options.phase), '-v7.3');
%
%     options.phase = 'EO';
%     fileName = makePSDFileName('labels', 'AllSubjs', options);
%     save(fullfile(N.searchFolder_4arch_rej_ICcats, fileName), sprintf('TrialsInfo_%s_all',options.phase));
%
%     fileName = makePSDFileName('psd', 'AllSubjs', options);
%     save(fullfile(N.searchFolder_4arch_rej_ICcats, folderName, fileName), sprintf('PSD_%s_all_norm', options.phase), '-v7.3');
% else
%     error('Event not suitable for this analysis')
% end



