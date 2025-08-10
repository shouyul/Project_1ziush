function [data_path, fig_path] = makeClassifierArchitecture(data_root, cfg)

switch cfg.class.condAnalysis
    case 'separate'
        % Tumbler
        folder_lvl1 = 'ConditionWise';
        folder_lvl2 = sprintf('%sPrediction', cfg.class.contrast);
        folder_lvl3 = '';
    case ''
        % VEP
        folder_lvl1 = '';
        folder_lvl2 = 'TrialTypePrediction';
        
        if strcmp(cfg.class.contrast, 'TrialType')
            folder_lvl3 = 'SixWay';
        elseif contains(cfg.class.contrast, 'Pairwise')
            folder_lvl3 = 'Pairwise';
        else
            folder_lvl3 = cfg.class.contrast;
        end
    otherwise
        error('Not coded yet')
end

folder_lvl4 = cfg.feat.folderName;
folder_lvl5 = cfg.class.model;

switch cfg.feat.type
    case 'psd'
        opts = cfg.psd;
        opts.SR = cfg.resample_freq;
        folder_lvl0 = makePSDFolderName(opts);
        arch = folder_lvl0;
        if ~isempty(folder_lvl1)
            arch = fullfile(arch,folder_lvl1);
        end
        arch = fullfile(arch,folder_lvl2);
        if ~isempty(folder_lvl3)
            arch = fullfile(arch,folder_lvl3);
        end
        arch = fullfile(arch,folder_lvl4);
        arch = fullfile(arch,folder_lvl5);
        %         if isempty(folder_lvl1)
        %             arch = fullfile(folder_lvl0, folder_lvl2, folder_lvl4, folder_lvl5);
        %         else
        %             arch = fullfile(folder_lvl0, folder_lvl1, folder_lvl2, folder_lvl4, folder_lvl5);
        %         end
    case 'amp'
        if ~isempty(folder_lvl1)
            arch = folder_lvl1;
        else
            arch = fullfile(folder_lvl1,folder_lvl2);
        end
        
        if ~isempty(folder_lvl3)
            arch = fullfile(arch,folder_lvl3);
        end
        arch = fullfile(arch,folder_lvl4);
        arch = fullfile(arch,folder_lvl5);
        %         if isempty(folder_lvl1)
        %             arch = fullfile(folder_lvl0, folder_lvl2, folder_lvl4, folder_lvl5);
        %         else
        %             arch = fullfile(folder_lvl0, folder_lvl1, folder_lvl2, folder_lvl4, folder_lvl5);
        %         end
        %         if isempty(folder_lvl1)
        %             arch = fullfile(folder_lvl2, folder_lvl4, folder_lvl5);
        %         else
        %             arch = fullfile(folder_lvl1, folder_lvl2, folder_lvl4, folder_lvl5);
        %         end
end

data_path = fullfile(data_root, arch);
if ~exist(data_path, 'dir')
    mkdir(data_path);
end
data_path = [data_path, filesep];

temp = strfind(data_root,'-subject-analysis');
arch2 = data_root(temp+length('-subject-analysis')+1:end);
fig_path = fullfile(cfg.figures_folder, 'Classifiers', arch2, arch);
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end
fig_path = [fig_path, filesep];
end