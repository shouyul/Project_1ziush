function EEG = events_check_VEP(EEG, cfg)
% Checks that all expected events are present for each trial
% 
% Output: Creates EEG.etc.TrialsInspection

subject = cfg.subjects(cfg.current_subject).id;
%% Specifics of the experiment
if cfg.subjects(cfg.current_subject).date < datetime(2023,7,30)
    blockLength = 12;
else
    blockLength = 36;
end
n_Blocks = length(unique([EEG.event.block]));
n_Trials = n_Blocks*blockLength;

%% Events
events = EEG.event;
evts_noBounds = events(~(strcmp({events.type},'boundary') | strcmp({events.type},'ExpEnd')));
%missing_base_start = false;
% struct markers with the BaseStart

%% Prepare columns of the dataset
%Condition = cell(n_Trials,1);
TrialInd = nan(n_Trials,1);
BlockInd = nan(n_Trials,1);
TrialType = cell(n_Trials,1);
% Latencies
Trial_urevent_seq = nan(n_Trials,2);
% EEG Completeness
Trial_perc = nan(n_Trials,1);
MaxBufferPre = nan(n_Trials,1);
MaxBufferPost = nan(n_Trials,1);
% ClosedEyes_perc = nan(n_Trials,1);
% Observation_perc = nan(n_Trials,1);
% Answer_perc = nan(n_Trials,1);

%% Create useful variables
startTrialEvents = strcmp({evts_noBounds.type},'TrialStart');
% obsStartEvents = strcmp({evts_noBounds.type},'OpenEyesSoundOn');
% obsEndEvents = strcmp({evts_noBounds.type},'ClosedEyesSoundOn');
% questionStartEvents = strcmp({evts_noBounds.type},'QuestionSoundOn');
% answerEvents = strcmp({evts_noBounds.type},'KeyPressed');
endTrialEvents = strcmp({evts_noBounds.type},'TrialStop');

%% Loop over trials
for tr = 1:n_Trials
    % useful variables
    bl = ceil(tr/blockLength);
    tr_inBl = tr - (bl-1)*blockLength;
    trialEvents = [evts_noBounds.block] == bl & [evts_noBounds.trial] == tr_inBl;
    
    % Fill columns
    TrialInd(tr) = tr_inBl;
    BlockInd(tr) = bl;
    
    % Check whether we are missing the trial
    if ~any(trialEvents & startTrialEvents)
        warning('Could not find trial start for Trial %d in Block %d.', tr_inBl, bl);
        if any(trialEvents)
            % Some events exist for this trial
            disp('Consider choosing a replacement for the start event...');
            % To do
            %disp('...');
        else
            disp('This trial was not recorded...');
            disp('Cannot infer trial type...')
            continue
        end
    else
        trialStartEvent = trialEvents & startTrialEvents;
        TrialType{tr} = evts_noBounds(trialStartEvent).stimulus;
        % Fill urevent sequence
        Trial_urevent_seq(tr,1) = evts_noBounds(trialStartEvent).urevent;
    end
    %Condition{tr} = evts_noBounds(trialStartEvent).condition;
    
    
    
%     if any(trialEvents & obsStartEvents)
%         Trial_urevent_seq(tr,2) = evts_noBounds(trialEvents & obsStartEvents).urevent;
%     else
%         warning('Could not find obs start for Trial %d in Block %d.', tr_inBl, bl);
%     end
%     
%     if any(trialEvents & obsEndEvents)
%         Trial_urevent_seq(tr,3) = evts_noBounds(trialEvents & obsEndEvents).urevent;
%     else
%         warning('Could not find obs end for Trial %d in Block %d.', tr_inBl, bl);
%     end
%     
%     if any(trialEvents & questionStartEvents)
%         Trial_urevent_seq(tr,4) = evts_noBounds(trialEvents & questionStartEvents).urevent;
%     else
%         warning('Could not find question start for Trial %d in Block %d.', tr_inBl, bl);
%     end
%     
%     if any(trialEvents & answerEvents)
%         if sum(trialEvents & answerEvents) > 1
%             warning('More than 1 answer event for Trial %d in Block %d.', tr_inBl, bl);
%         else
%             Trial_urevent_seq(tr,5) = evts_noBounds(trialEvents & answerEvents).urevent;
%         end
%     else
%         warning('Could not find answer for Trial %d in Block %d.', tr_inBl, bl);
%     end
    
    if any(trialEvents & endTrialEvents)
        if ~isnan(Trial_urevent_seq(tr,1)) % not bothering including an end event is start event wasn't found
            Trial_urevent_seq(tr,2) = evts_noBounds(trialEvents & endTrialEvents).urevent;
        end
    else
        warning('Could not find trial end for Trial %d in Block %d.', tr_inBl, bl);
    end
end

%% Second loop for EEG completeness
for tr = 1:n_Trials
    bl = ceil(tr/blockLength);
    tr_inBl = tr - (bl-1)*blockLength;
    
    % Full trial first
    limits = [evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,1)).latency,...
        evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,end)).latency];
    
    if isempty(limits)
        warning('Not enough latency information for trial %d', tr);
        Trial_perc(tr) = 0;
    elseif limits(1) >= limits(2)
        error('Wrong order of events for trial %d', tr);
    else        
        if isnan(Trial_perc(tr))
            Trial_perc(tr) = 100 - computeMissingEEGData(EEG,limits(1),limits(2));
            % Otherwise this has already been set by a special case
        end
        
%         if Trial_perc(tr) == 100
%             ClosedEyes_perc(tr) = 100;
%             Observation_perc(tr) = 100;
%             Answer_perc(tr) = 100;
%         elseif Trial_perc(tr) > 0
%             
%             limits_EC = [evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,1)).latency,...
%                 evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,2)).latency];
%             if isempty(limits_EC)
%                 warning('Not enough latency information for trial %d, in Closed Eyes Baseline', tr);
%                 ClosedEyes_perc(tr) = 0;
%             elseif limits_EC(1) >= limits_EC(2)
%                 error('Wrong order of events for trial %d, in Closed Eyes Baseline', tr);
%             else
%                 ClosedEyes_perc(tr) = 100 - computeMissingEEGData(EEG, limits_EC(1),limits_EC(2));
%             end
%             
%             limits_Obs = [evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,2)).latency,...
%                 evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,3)).latency];
%             if isempty(limits_Obs)
%                 warning('Not enough latency information for trial %d, in Observation period', tr);
%                 Observation_perc(tr) = 0;
%             elseif limits_Obs(1) >= limits_Obs(2)
%                 error('Wrong order of events for trial %d, in Observation period', tr);
%             else
%                 Observation_perc(tr) = 100 - computeMissingEEGData(EEG, limits_Obs(1),limits_Obs(2));
%             end
%             
%             limits_Ans = [evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,3)).latency,...
%                 evts_noBounds([evts_noBounds.urevent] == Trial_urevent_seq(tr,end)).latency];
%             if isempty(limits_Ans)
%                 warning('Not enough latency information for trial %d, in Answer period', tr);
%                 Answer_perc(tr) = 0;
%             elseif limits_Ans(1) >= limits_Ans(2)
%                 error('Wrong order of events for trial %d, in Answer period', tr);
%             else
%                 Answer_perc(tr) = 100 - computeMissingEEGData(EEG, limits_Ans(1),limits_Ans(2));
%             end
%             
%         else
%             ClosedEyes_perc(tr) = 0;
%             Observation_perc(tr) = 0;
%             Answer_perc(tr) = 0;
%         end
    end
    
    if Trial_perc(tr) == 100
        start_ev = find([EEG.event.urevent] == Trial_urevent_seq(tr,1));
        MaxBufferPre(tr) = (EEG.times(limits(1)) - EEG.times(searchNextNaN(EEG, start_ev, 'backward')))/1000;
        stop_ev = find([EEG.event.urevent] == Trial_urevent_seq(tr,end));
        MaxBufferPost(tr) = (EEG.times(searchNextNaN(EEG, stop_ev, 'forward')) - EEG.times(limits(2)))/1000;
    end
end

EEGCompletenessSummary = table(BlockInd, TrialInd, TrialType, Trial_urevent_seq, Trial_perc, MaxBufferPre, MaxBufferPost);

%% Plot Summary info:
if ~exist(fullfile(cfg.figures_folder, 'MissingData'),'dir')
    mkdir(fullfile(cfg.figures_folder, 'MissingData'));
end

figure
for b = 1:n_Blocks
    subplot(ceil(n_Blocks/2), 2, b)
    hold on
    xlim([0,blockLength+1])
    trials_block = BlockInd == b;
%     trials_control_bl = trials_block & strcmp(TrialType,'Control');
%     trials_disc_bl = trials_block & strcmp(TrialType,'Disc');
%     trials_3rings_bl = trials_block & strcmp(TrialType,'3Rings');
%     trials_5rings_bl = trials_block & strcmp(TrialType,'5Rings');
%     trials_45bar_bl = trials_block & strcmp(TrialType,'45Bar');
%     trials_min45bar_bl = trials_block & strcmp(TrialType,'-45Bar');
    
    bar(1:blockLength, Trial_perc(trials_block));
%     bar(TrialInd(trials_control_bl)', Trial_perc(trials_control_bl)',0.5);
%     bar(TrialInd(trials_disc_bl)', Trial_perc(trials_disc_bl)',0.5);
%     bar(TrialInd(trials_3rings_bl)', Trial_perc(trials_3rings_bl)',0.5);
%     bar(TrialInd(trials_5rings_bl)', Trial_perc(trials_5rings_bl)',0.5);
%     bar(TrialInd(trials_45bar_bl)', Trial_perc(trials_45bar_bl)',0.5);
%     bar(TrialInd(trials_min45bar_bl)', Trial_perc(trials_min45bar_bl)',0.5);
    title(sprintf('Block %d', b))
    ylim([0,100])
    xlabel('Trial')
    ylabel('% Complete')
%     if b == 1
%         legend({'Control', 'Disc', '3Rings', '5Rings', '45Bar', '-45Bar'}, 'Position', [0.825,0.925,0.1,0.05])
%     end
end
subject = cfg.subjects(cfg.current_subject).id;
suptitle(subject);
saveCurrentFig([cfg.figures_folder 'MissingData' filesep],...
    [subject, '_TrialsInspectionResults'], {'png'}, [1400,700]);

% figure
% for b = 1:n_Blocks
%     subplot(n_Blocks/2, 2, b)
%     trials_bl = BlockInd == b;
%     
%     bar(TrialInd(trials_bl), table2array(EEGCompletenessSummary(trials_bl,9:11)))
%     xticks(TrialInd(trials_bl));
%     xticklabels(TrialType(trials_bl))
%     xtickangle(45)
%     title(sprintf('Block %d', b))
%     ylim([0,100])
%     ylabel('% Complete')
%     if b == 1
%         legend({'Eyes Closed', 'Observations', 'Answer'}, 'Position', [0.825,0.93,0.1,0.05])
%     end
% end
% subject = cfg.subjects(cfg.current_subject).id;
% suptitle(subject);
% saveCurrentFig([cfg.figures_folder 'MissingData' filesep],...
%     [subject, '_TrialsInspectionResults_Details'], {'png'}, [1400,700]);

EEG.etc.TrialsInspection = EEGCompletenessSummary;

%% Helper Functions ------------------------------------------
    function  perc = computeMissingEEGData(EEG,lat_start,lat_stop)
        pnts = lat_start:lat_stop;
        EEG_chans = strcmp({EEG.chanlocs.type}, 'EEG');
        data2inspect = EEG.data(EEG_chans,pnts);
        NanInspection = isnan(data2inspect);
        NanInspection = sum(NanInspection,1)>0;
        perc = 100*sum(NanInspection)/length(NanInspection);
    end

    function lat = searchNextNaN(EEG, ref_ind, direction)
        switch direction
            case 'backward'
                data2inspect = EEG.data(:,1:(EEG.event(ref_ind).latency-1));
            case 'forward'
                data2inspect = EEG.data(:,(EEG.event(ref_ind).latency+1):end);
        end
        
        NanInspection = isnan(data2inspect);
        NanInspection = sum(NanInspection,1)>0;
        
        switch direction
            case 'backward'
                next_bound = find(strcmp({EEG.event(1:(ref_ind-1)).type}, 'boundary'),1,'last');
                if isempty(next_bound)
                    lat = find(NanInspection, 1,'last');
                    if isempty(lat)
                        % No NaNs found before the event latency, probably
                        % the start of the recording (with no boundary
                        % event starting the recording)
                        lat = 1;
                    end
                else
                    lat = max(EEG.event(next_bound).latency, find(NanInspection,1,'last'));
                end
            case 'forward'
                next_bound = find(strcmp({EEG.event((ref_ind+1):end).type}, 'boundary'),1);
                if isempty(next_bound)
                    lat = EEG.event(ref_ind).latency + find(NanInspection,1);
                    if isempty(lat)
                        % No NaNs found after the event latency, probably
                        % the end of the recording (with no boundary
                        % event ending the recording)
                        lat = EEG.pnts;
                    end
                else
                    lat = min(EEG.event(ref_ind+next_bound).latency, EEG.event(ref_ind).latency + find(NanInspection,1));
                end
        end
        % Make sure the output is an integer
        lat = round(lat);
    end
end