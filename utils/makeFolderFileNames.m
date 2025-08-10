function names = makeFolderFileNames(cfg, subject)
names.searchFolder_1 = [cfg.study_folder cfg.raw_EEGLAB_data_folder subject filesep];
createFolder(names.searchFolder_1);
if isempty(cfg.channel_selection)
    names.searchFolder_2 = [cfg.study_folder cfg.preprocessing_folder];
else
    names.searchFolder_2 = [cfg.study_folder cfg.preprocessing_folder 'subSelect' filesep];
end
createFolder(names.searchFolder_2);
names.searchFolder_2arch = [names.searchFolder_2 lower(cfg.globalArchitecture) filesep];
createFolder(names.searchFolder_2arch);
if isempty(cfg.channel_selection)
    names.searchFolder_3 = [cfg.study_folder cfg.single_subject_analysis_folder];
else
    names.searchFolder_3 = [cfg.study_folder cfg.single_subject_analysis_folder 'subSelect' filesep];
end
createFolder(names.searchFolder_3);
names.searchFolder_3arch = [names.searchFolder_3 lower(cfg.globalArchitecture) filesep];
createFolder(names.searchFolder_3arch);
if isempty(cfg.channel_selection)
    names.searchFolder_4 = [cfg.study_folder cfg.multi_subject_analysis_folder];
else
    names.searchFolder_4 = [cfg.study_folder cfg.multi_subject_analysis_folder 'subSelect' filesep];
end
createFolder(names.searchFolder_4);
names.searchFolder_4arch = [names.searchFolder_4 lower(cfg.globalArchitecture) filesep];
createFolder(names.searchFolder_4arch);
switch cfg.badSampsRejection
    case 'manual'
        names.searchFolder_2arch_rej = [names.searchFolder_2arch 'Manual' filesep];
        names.searchFolder_3arch_rej = [names.searchFolder_3arch 'Manual' filesep];
        names.searchFolder_4arch_rej = [names.searchFolder_4arch 'Manual' filesep];
    case 'app'
        names.searchFolder_2arch_rej = [names.searchFolder_2arch 'APP' filesep];
        names.searchFolder_3arch_rej = [names.searchFolder_3arch 'APP' filesep];
        names.searchFolder_4arch_rej = [names.searchFolder_4arch 'APP' filesep];
    case 'asr'
        switch cfg.ASR_use
            case 'rewrite'
                names.searchFolder_2arch_rej = [names.searchFolder_2arch 'ASR_corrected' filesep];
                names.searchFolder_3arch_rej = [names.searchFolder_3arch 'ASR_corrected' filesep];
                names.searchFolder_4arch_rej = [names.searchFolder_4arch 'ASR_corrected' filesep];
            case 'reject'
                names.searchFolder_2arch_rej = [names.searchFolder_2arch 'ASR_rejected' filesep];
                names.searchFolder_3arch_rej = [names.searchFolder_3arch 'ASR_rejected' filesep];
                names.searchFolder_4arch_rej = [names.searchFolder_4arch 'ASR_rejected' filesep];
            otherwise
                error('Unknown ASR use')
        end
    otherwise
        error('Unknown bad samples rejection')
end
createFolder(names.searchFolder_2arch_rej);
createFolder(names.searchFolder_3arch_rej);
createFolder(names.searchFolder_4arch_rej);

if strcmp(cfg.ICselection, 'manual')
    for c = 1:numel(cfg.cats2keep)
        if c == 1
            names.searchFolder_2arch_rej_ICcats = [names.searchFolder_2arch_rej cfg.cats2keep{c}];
            names.searchFolder_3arch_rej_ICcats = [names.searchFolder_3arch_rej cfg.cats2keep{c}];
            names.searchFolder_4arch_rej_ICcats = [names.searchFolder_4arch_rej cfg.cats2keep{c}];
        else
            names.searchFolder_2arch_rej_ICcats = [names.searchFolder_2arch_rej_ICcats '_' cfg.cats2keep{c}];
            names.searchFolder_3arch_rej_ICcats = [names.searchFolder_3arch_rej_ICcats '_' cfg.cats2keep{c}];
            names.searchFolder_4arch_rej_ICcats = [names.searchFolder_4arch_rej_ICcats '_' cfg.cats2keep{c}];
        end
    end
else
    names.searchFolder_2arch_rej_ICcats = [names.searchFolder_2arch_rej 'BrainAuto'];
    names.searchFolder_3arch_rej_ICcats = [names.searchFolder_3arch_rej 'BrainAuto'];
    names.searchFolder_4arch_rej_ICcats = [names.searchFolder_4arch_rej 'BrainAuto'];
end
names.searchFolder_2arch_rej_ICcats = [names.searchFolder_2arch_rej_ICcats filesep];
createFolder(names.searchFolder_2arch_rej_ICcats);
names.searchFolder_3arch_rej_ICcats = [names.searchFolder_3arch_rej_ICcats filesep];
createFolder(names.searchFolder_3arch_rej_ICcats);
names.searchFolder_4arch_rej_ICcats = [names.searchFolder_4arch_rej_ICcats filesep];
createFolder(names.searchFolder_4arch_rej_ICcats);

names.postimportFile = [subject '_EEG_' cfg.merged_filename];
names.postimport_withMoCapFile = [subject '_AllStreams_' cfg.merged_filename];
names.preparedFile = [subject '_' cfg.prepared_filename];
names.prepared_withMoCapFile = [subject '_AllStreams_' cfg.prepared_filename];
names.MoCapPreparedFile = [subject '_MoCapPrepared.set'];
names.nobadchansFile = [subject '_' cfg.BadChansRemoved_filename];
names.preICAFile = [subject '_' cfg.beforeICA_filename];
names.postICAFile = [subject '_' cfg.icaOutput_filename];
names.dipfitFile = [subject '_' cfg.dipolesFitted_filename];
names.IClabelledFile = [subject '_' cfg.icLabelled_filename];
names.postLabelingFile = [subject '_' cfg.icaSelect_filename];
names.finalFile = [subject '_' cfg.postICA_tempRej_filename];

switch cfg.epochs.event
    case 'EyesOpening'
        prefix_epoch = sprintf('%s_eyeOp',subject);
    case 'TumblerVisible'
        prefix_epoch = sprintf('%s_TumbVis',subject);
    case 'TrialStart'
        prefix_epoch = sprintf('%s_TrStart',subject);
    otherwise
        error('Unknown event type')
end

switch lower(cfg.epochs.window)
    case 'full'
        names.epochedFile = sprintf('%s_full_%s',prefix_epoch,cfg.epoched_filename);
    case 'fixed'
        if cfg.epochs.limits_wdw(1) <= 0
            if all(round(cfg.epochs.limits_wdw)==cfg.epochs.limits_wdw)
                names.epochedFile = sprintf('%s_bef%d-aft%d_%s',prefix_epoch,...
                    abs(cfg.epochs.limits_wdw(1)), cfg.epochs.limits_wdw(2), cfg.epoched_filename);
            else
                names.epochedFile = sprintf('%s_bef%dms-aft%dms_%s',prefix_epoch,...
                    abs(cfg.epochs.limits_wdw(1))*1000, cfg.epochs.limits_wdw(2)*1000, cfg.epoched_filename);
            end
        else
            if all(round(cfg.epochs.limits_wdw)==cfg.epochs.limits_wdw)
                names.epochedFile = sprintf('%s_aft%d-aft%d_%s',prefix_epoch,...
                    abs(cfg.epochs.limits_wdw(1)), cfg.epochs.limits_wdw(2), cfg.epoched_filename);
            else
                names.epochedFile = sprintf('%s_aft%dms-aft%dms_%s',prefix_epoch,...
                    abs(cfg.epochs.limits_wdw(1))*1000, cfg.epochs.limits_wdw(2)*1000, cfg.epoched_filename);
            end
        end
    otherwise
        error('Unknown window type')
end
names.baselineEpochedFile = sprintf('%s_%s',subject,cfg.base_epoched_filename);

    function createFolder(path)
        if ~exist(path, 'dir')
            mkdir(path);
        end
    end
end