function [main, raw, elec, raw_EEGLAB, preproc, SS_ana, MS_ana, Figs, MetaFile, MetaTab] =...
    getMainFoldersNames_PIONEER(user, task)
% Define here the folders access.
% Inputs:
%   user        - permits to specify the user to adapt the datapath
%   task        - permits to specify the task to adapt the datapath

switch lower(user)
    case 'sl'
        main = fullfile('D:\Data_PIONEER',sprintf('%sTask',task),'analysis', filesep);
        Figs = fullfile('D:\Data_PIONEER',sprintf('%sTask',task),'figures', filesep);
        MetaFile = fullfile('D:\Data_PIONEER', sprintf('%sTask',task),...
            sprintf('TrackingEEGanalysis_PIONEER_%s.xlsx',task));
        MetaTab = 'SubjectsInfo';
    case 'jb'
        main = fullfile('/Users/jean-baptiste/Documents/DATA/PIONEER-EEG/ANALYSIS',sprintf('%sTask',task),'analysis', filesep);
        Figs = fullfile('/Users/jean-baptiste/Documents/DATA/PIONEER-EEG/ANALYSIS',sprintf('%sTask',task),'figures', filesep);
        MetaFile = fullfile('/Users/jean-baptiste/Documents/DATA/PIONEER-EEG/ANALYSIS', sprintf('%sTask',task),...
            sprintf('TrackingEEGanalysis_PIONEER_%s.xlsx',task));
        MetaTab = 'SubjectsInfo';
    otherwise
        error('Unknown user');
end
raw = ['0_raw-data' filesep];
elec = ['0_electrodes' filesep]; % Useful if you have electrodes digitization data
raw_EEGLAB = ['1_raw-eeglab' filesep];
preproc = ['2_preprocessing' filesep];
SS_ana = ['3_single-subject-analysis' filesep];
MS_ana = ['4_multi-subject-analysis' filesep];
end

