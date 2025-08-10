function EEG = BehaviorReportTumbler(EEG, cfg)
% Should be called after fillNaNs and events_check but before resampling or selection
% if all the trials are present
%
% Behavioral reports specific to the experiment
%
% Output: Creates EEG.etc.TrialData and also saves it externally

events = EEG.event;
evts_noBounds = events(~strcmp({events.type},'boundary'));
answerEvents = strcmp({evts_noBounds.type},'KeyPressed');
confRateEvents = strcmp({evts_noBounds.type},'ConfidenceRating');

DataBase = struct('ID', [], 'Gender', [], 'Block', [], 'Condition', [], 'Trial', [], 'TrialType', [],...
    'Answer', [], 'Delay_ms', [], 'ConfidenceLevel', [], 'urevent_seq', [], 'EEGComplete_trial',[], 'EEGComplete_details',[]);

% Subject specific variables
subject = cfg.subjects(cfg.current_subject).id;
gender = cfg.subjects(cfg.current_subject).gender;

TrialInspection = EEG.etc.TrialsInspection;

i=1;
for tr = 1:size(TrialInspection,1)
    bl = TrialInspection.BlockInd(tr);
    tr_inBl = TrialInspection.TrialInd(tr);
    trialEvents = [evts_noBounds.block] == bl & [evts_noBounds.trial] == tr_inBl;
    
    if strcmp(subject, 'P1001-3') && tr == 40
        % Special case from experimental feedback
        asw = 'Present';
        EEG.event(...
            findInStructWithEmpties(EEG.event,'urevent',evts_noBounds(trialEvents & answerEvents).urevent)...
            ).answer = asw;
    else
        asw = evts_noBounds(trialEvents & answerEvents).answer;
    end
    
    if any(trialEvents & confRateEvents)        
        cfd_lvl = str2double(evts_noBounds(trialEvents & confRateEvents).answer);
    else
        cfd_lvl = NaN;
    end
    
    % Concatenate EEG completness data
    eeg_perc = [TrialInspection.ClosedEyes_perc(tr),...
        TrialInspection.Observation_perc(tr),...
        TrialInspection.Answer_perc(tr)];
    
    %fprintf('Block %d, Trial %d.\n', bl, tr_inBl)
    DataBase(i) = struct('ID', subject, 'Gender', gender,...
        'Block', bl, 'Condition', TrialInspection.Condition(tr),...
        'Trial', tr_inBl, 'TrialType', TrialInspection.TrialType(tr),...
        'Answer', asw,...
        'Delay_ms', evts_noBounds(trialEvents & answerEvents).delay*1000,...
        'ConfidenceLevel', cfd_lvl,...
        'urevent_seq', TrialInspection.Trial_urevent_seq(tr,:),...
        'EEGComplete_trial', TrialInspection.Trial_perc(tr),...
        'EEGComplete_details', eeg_perc);
    i=i+1;
end

DataBase = struct2table(DataBase);

%% Plot Behavioral data (confusion matrix)
if ~exist(fullfile(cfg.figures_folder, 'Behavior'), 'dir')
    mkdir(fullfile(cfg.figures_folder, 'Behavior'));
end

%Per condition
for c = 1:2
    conf_mat = zeros(2,2);
    if c == 1
        DataCond = DataBase(contains(DataBase.Condition, 'ON'),:);
        titleCond = 'GogglesON';
    else
        DataCond = DataBase(contains(DataBase.Condition, 'OFF'),:);
        titleCond = 'GogglesOFF';
    end
    
    for tr = 1:size(DataCond,1)
        switch DataCond.TrialType{tr}
            case 'WithObject'
                col = 1;
            case 'WithoutObject'
                col = 2;
        end
        
        switch DataCond.Answer{tr}
            case 'Present'
                ln = 1;
            case 'Absent'
                ln = 2;
        end
        
        conf_mat(col,ln) = conf_mat(col,ln)+1;
    end
    
    figure;
    plotConfMat(conf_mat', {'Obj Present','Obj Absent'});
    suptitle(sprintf('%s - %s', subject, titleCond));
    saveCurrentFig([cfg.figures_folder 'Behavior' filesep],...
        sprintf('%s_ConfMat_%s', subject, titleCond), {'png'}, []);
end

%Per block
blocks = unique(DataBase.Block);
for b = blocks'
    conf_mat = zeros(2,2);
    DataBlock = DataBase(DataBase.Block == b,:);
    titleBlock = sprintf('Block%d',b);
    blockType = DataBlock.Condition{1};
    
    for tr = 1:size(DataBlock,1)
        switch DataBlock.TrialType{tr}
            case 'WithObject'
                col = 1;
            case 'WithoutObject'
                col = 2;
        end
        
        switch DataBlock.Answer{tr}
            case 'Present'
                ln = 1;
            case 'Absent'
                ln = 2;
        end
        
        conf_mat(col,ln) = conf_mat(col,ln)+1;
    end
    
    figure;
    plotConfMat(conf_mat', {'Obj Present','Obj Absent'});
    suptitle(sprintf('%s - %s (%s)', subject, titleBlock, blockType));
    saveCurrentFig([cfg.figures_folder 'Behavior' filesep],...
        sprintf('%s_ConfMat_%s', subject, titleBlock), {'png'}, []);
end

%% Save the output structure
EEG.etc.TrialData = DataBase;
dirname = [cfg.study_folder cfg.preprocessing_folder];
fname = [subject '_DataBase.mat'];
save(fullfile(dirname, fname), 'DataBase')
end