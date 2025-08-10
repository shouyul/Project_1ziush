
config_PIONEER_Tumbler
subject = study_config.subjects(subject_ind).id;
cprintf('*Blue',[subject '\n']);

load([study_config.study_folder study_config.preprocessing_folder subject '_DataBase.mat']);
%load('D:\Data_PIONEER\TumblerTask\analysis\2_preprocessing\P1002-3_DataBase.mat')
trials_ggON = contains(DataBase.Condition,'ON');
trials_WO = strcmp(DataBase.TrialType,'WithObject');
trials_ansPre = strcmp(DataBase.Answer,'Present');

%figure
% violin([DataBase.ConfidenceLevel(trials_ggON & trials_WO),...
%     DataBase.ConfidenceLevel(trials_ggON & ~trials_WO)],...
%     'facecolor', [0,1,0;1,0,0], 'facealpha', 0.25,...
%     'xlabel', {'WithObject','WithoutObject'}, 'bw', 0.5);
% hold on
% scatter(0.95+0.1*rand(sum(trials_ggON & trials_WO),1),...
%     DataBase.ConfidenceLevel(trials_ggON & trials_WO),32,[0,1,0]);
% scatter(1.95+0.1*rand(sum(trials_ggON & ~trials_WO),1),...
%     DataBase.ConfidenceLevel(trials_ggON & ~trials_WO),32,[1,0,0]);
% ylabel('Confidence Rating')

figure
% violin({DataBase.ConfidenceLevel(trials_ggON & trials_ansPre),...
%     DataBase.ConfidenceLevel(trials_ggON & ~trials_ansPre)},...
%     'facecolor', [0,0,1;1,1,1], 'facealpha', 0.25,...
%     'xlabel', {'Present','Absent'}, 'bw', 0.5);
% hold on
% scatter(0.95+0.1*rand(sum(trials_ggON & trials_ansPre),1),...
%     DataBase.ConfidenceLevel(trials_ggON & trials_ansPre),32,[0,0,1]);
% scatter(1.95+0.1*rand(sum(trials_ggON & ~trials_ansPre),1),...
%     DataBase.ConfidenceLevel(trials_ggON & ~trials_ansPre),32,[0,0,0]);
% ylabel('Confidence Rating')

edges = [0.5:1:5.5];
N_pre = histcounts(DataBase.ConfidenceLevel(trials_ggON & trials_ansPre),edges);
N_abs = histcounts(DataBase.ConfidenceLevel(trials_ggON & ~trials_ansPre),edges);
bar([N_pre;N_abs]')
ylim([0,15])
xlabel('Confidence Rating')
ylabel('Number of trials')
legend({sprintf('Present (N=%d)',sum(N_pre)), sprintf('Absent (N=%d)',sum(N_abs))})
title('Confidence rating per answer type')

saveCurrentFig([study_config.figures_folder 'Behavior' filesep], sprintf('%s_Confid_by_Answer-GogglesON', subject), {'png','svg'}, []);


figure
% violin({DataBase.ConfidenceLevel(trials_ggON & trials_WO & trials_ansPre),...
%     DataBase.ConfidenceLevel(trials_ggON & ~trials_WO & trials_ansPre)},...
%     'facecolor', [0,1,1;1,0,1], 'facealpha', 0.25,...
%     'xlabel', {'Hit','FalsePos'}, 'bw', 0.5);
% hold on
% scatter(0.95+0.1*rand(sum(trials_ggON & trials_WO & trials_ansPre),1),...
%     DataBase.ConfidenceLevel(trials_ggON & trials_WO & trials_ansPre),32,[0,1,1]);
% scatter(1.95+0.1*rand(sum(trials_ggON & ~trials_WO & trials_ansPre),1),...
%     DataBase.ConfidenceLevel(trials_ggON & ~trials_WO & trials_ansPre),32,[1,0,1]);
% ylabel('Confidence Rating')

N_hits = histcounts(DataBase.ConfidenceLevel(trials_ggON & trials_WO & trials_ansPre),edges);
N_FA = histcounts(DataBase.ConfidenceLevel(trials_ggON & ~trials_WO & trials_ansPre),edges);
bar([N_hits;N_FA]')
ylim([0,15])
xlabel('Confidence Rating')
ylabel('Number of trials')
legend({sprintf('Hits (N=%d)',sum(N_hits)), sprintf('False Alarms (N=%d)',sum(N_FA))})
title('Confidence rating per correct response when answered Present')

saveCurrentFig([study_config.figures_folder 'Behavior' filesep], sprintf('%s_Confid_by_Correctness-ObjPresent-GogglesON', subject),...
    {'png','svg'}, []);



