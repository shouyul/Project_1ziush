function EEG = preprocess(EEG, config, doBadChans)
% doBadChans: specific to the 'bemobil' pipeline. Avoid recomputing bad
% channels. This takes into account that the EEG dataset has already
% undergone this step. The bad channels corrected dataset will not be
% loaded within this function.

subject = config.subjects(config.current_subject).id;
N = makeFolderFileNames(config, subject);

switch lower(config.globalArchitecture)
    case 'simple'
        error('This pipeline is not up to date check how to retrieve bad electrodes.')
        % Remove bad channels flagged during recording
        if isfield(config.badElectrodes, subject)
            badElectrodes = config.badElectrodes.(subject);
            bads = [];
            for ch = 1:numel(badElectrodes)
                bads = [bads, find(contains({noisyOut.channelLocations.labels}, badElectrodes{ch}))];
            end
            EEG = pop_select(EEG,'nochannel', bads);
        end
        
        % Filtering between freq limits
        EEG_filt = custom_filter(EEG, config.low_cut_off, config.high_cut_off);
        % To enable later vizualization of ASR effects:
        pop_saveset(EEG_filt, 'filename', [subject '_' config.ASRin_filename],...
            'filepath', [config.study_folder config.preprocessing_folder ...
            lower(config.globalArchitecture) filesep 'ASR']);
        
        %%%%%%%%%% clean_artifacts: All-in-one function for artifact removal, including ASR.
        % This function removes flatline channels, low-frequency drifts, noisy channels, short-time bursts
        % and incompletely repaired segments from the data.
        % Tip: Any of the core parameters can also be passed in as [] to use the respective default of the underlying functions, or as 'off' to disable
        % it entirely.
        EEG_filt.etc.ASRBurstCrit = config.burst_crit;
        switch config.ASR_use
            case 'reject'
                % Use the function for channel rejection and very bad temporal segments only
                [EEG_asr,~,~] = clean_artifacts(EEG_filt, 'Highpass', 'off',...
                    'BurstCriterion', config.burst_crit, 'BurstRejection', 'off');
                
                % To enable later vizualization of ASR effects:
                pop_saveset(EEG_asr, 'filename', [subject '_' config.ASRout_filename],...
                    'filepath', [config.study_folder config.preprocessing_folder ...
                    lower(config.globalArchitecture) filesep 'ASR_rejected' filesep 'ASR_output']);
                
                %Visualize the difference with:
                %vis_artifacts(EEG_asr,EEG_filt);
                
                if isfield(EEG_asr.etc,'clean_channel_mask')
                    EEG_preproc = pop_select(EEG_filt,'nochannel',find(EEG_asr.etc.clean_channel_mask==0));
                    EEG_preproc.etc.clean_channel_mask = EEG_asr.etc.clean_channel_mask;
                    
                    %Interpolate all the removed channels
                    EEG_preproc = pop_interp(EEG_preproc, EEG.chanlocs, 'spherical');
                else
                    EEG_preproc = EEG_filt;
                    EEG_preproc.etc.clean_channel_mask = true(EEG_preproc.nbchan,1);
                end
                
                if isfield(EEG_asr.etc,'clean_sample_mask') && sum(EEG_asr.etc.clean_sample_mask)<EEG_asr.pnts
                    sample_mask = EEG_asr.etc.clean_sample_mask;
                    if sample_mask(1)
                        interval_starts = [1, find(diff(sample_mask)==1)+1];
                    else
                        interval_starts = find(diff(sample_mask)==1)+1;
                    end
                    if sample_mask(end)
                        interval_ends = [find(diff(sample_mask)==-1),length(sample_mask)];
                    else
                        interval_ends = find(diff(sample_mask)==-1);
                    end
                    EEG_preproc = pop_select(EEG_preproc,'time',cat(2,interval_starts',interval_ends'));
                    EEG_preproc = eeg_checkset(EEG_preproc);
                    EEG_preproc.etc.clean_sample_mask = sample_mask;
                else
                    EEG_preproc.etc.clean_sample_mask = true(EEG_preproc.pnts,1);
                end
                
                clear EEG_asr
            case 'rewrite'
                [EEG_preproc,~,~] = clean_artifacts(EEG_filt, 'Highpass', 'off',...
                    'BurstCriterion', config.burst_crit, 'BurstRejection', 'off');
                
                % To enable later vizualization of ASR effects:
                pop_saveset(EEG_preproc, 'filename', [subject '_' config.ASRout_filename],...
                    'filepath', [config.study_folder config.preprocessing_folder ...
                    lower(config.globalArchitecture) filesep 'ASR_corrected' filesep 'ASR_output']);
                
                % Visualize the difference with:
                %vis_artifacts(EEG_preproc,EEG_filt);
                
                if isfield(EEG_preproc.etc,'clean_channel_mask')
                    %Interpolate all the removed channels
                    EEG_preproc = pop_interp(EEG_preproc, EEG.chanlocs, 'spherical');
                else
                    EEG_preproc.etc.clean_channel_mask = true(EEG_preproc.nbchan,1);
                end
        end
        
        % Re-reference the data to average
        EEG_preproc.nbchan = EEG_preproc.nbchan+1;
        EEG_preproc.data(end+1,:) = zeros(1, EEG_preproc.pnts);
        EEG_preproc.chanlocs(EEG_preproc.nbchan).labels = 'initialReference';
        EEG_preproc = pop_reref(EEG_preproc, []);
        EEG_preproc = pop_select(EEG_preproc,'nochannel',{'initialReference'});
        
        % Filtering included and propagated within this step
        EEG = EEG_preproc;
        
    case 'bemobil'
        if doBadChans
            % 1. Search for bad channels
            [noisyOut] = FindBadChannels(EEG, config);
            EEG.etc.noisyChannelsDetection = noisyOut;
            
            % 2. Interpolation:
            chans_to_interp = noisyOut.noisyChannels.all;            
            if ~isempty(chans_to_interp)
                disp('The following channels will be interpolated:')
                disp({EEG.chanlocs(chans_to_interp).labels})
                
                EOG_channels = find(strcmp({EEG.chanlocs.type},'EOG'));
                if ~isempty(EOG_channels)
                    % pop_interp doesn't exclude EOG channels from the
                    % good data when interpolating
                    % remove EOG channels
                    for ch = 1:length(chans_to_interp)
                        chans_to_interp(ch) = chans_to_interp(ch) - sum(EOG_channels < chans_to_interp(ch));
                    end                    
                    EEG_noEOG = pop_select(EEG, 'nochannel', EOG_channels);
                    % Interpolate
                    EEG_interp = pop_interp(EEG_noEOG, chans_to_interp, 'spherical');
                    % Replace EOG channels
                    EEG_interp.chanlocs = EEG.chanlocs;
                    EEG_interp.nbchan = EEG.nbchan;
                    datatmp = zeros(EEG.nbchan, EEG.pnts);
                    for ch = 1:EEG.nbchan
                        if sum(EOG_channels==ch)>0
                            datatmp(ch,:) = EEG.data(ch,:);
                        else
                            datatmp(ch,:) = EEG_interp.data(ch-sum(EOG_channels < ch),:);
                        end                         
                    end
                    EEG_interp.data = datatmp;            
                else
                    EEG_interp = pop_interp(EEG, chans_to_interp, 'spherical');
                end  
            else
                disp('No channels to interpolate.')
                EEG_interp = EEG;
            end
            clear EEG
            
            % 3. Rereferencing:
            EEG_interp.nbchan = EEG_interp.nbchan+1;
            EEG_interp.data(end+1,:) = zeros(1, EEG_interp.pnts);
            EEG_interp.chanlocs(EEG_interp.nbchan).labels = 'initialReference';
            EEG_interp_avRef = pop_reref(EEG_interp, []);
            EEG_interp_avRef = pop_select(EEG_interp_avRef,'nochannel',{'initialReference'});
            clear EEG_interp
            
            % Save this step
            pop_saveset(EEG_interp_avRef, 'filename', N.nobadchansFile, 'filepath', N.searchFolder_2arch);
        else
            % The EEG file should already be corrected for bad channels
            EEG_interp_avRef = EEG;
        end
        
        switch config.badSampsRejection
            case 'manual'
                EEG_manual = rejectBadTempsManually(EEG_interp_avRef, config);
                EEG = EEG_manual;
            case 'app'
                EEG_app = rejectBadTempsWithAPP(EEG_interp_avRef, config, false);
                EEG = EEG_app;
            case 'asr'
                switch config.ASR_use
                    case 'rewrite'
                        EEG_asr = rewriteBadTempsWithASR(EEG_interp_avRef, config);
                        EEG = EEG_asr;
                    case 'reject'
                        EEG_asr = rejectBadTempsWithASR(EEG_interp_avRef, config, false);
                        EEG = EEG_asr;
                end
            case 'autoMoBI'
                EEG_autoMoBI = rejectBadTempsWithAutoMoBI(EEG_interp_avRef, config);
                EEG = EEG_autoMoBI;
        end
end
end