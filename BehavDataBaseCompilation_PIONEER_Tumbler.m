clear all;
config_PIONEER_Tumbler;

%% Compile some general data for further analyses
% Options for this script
recompute = false;

%% Get the folders names right:
fig_path = fullfile(study_config.figures_folder, 'Epochs');
if ~exist(fig_path, 'dir')
    mkdir(fig_path)
end

N = makeFolderFileNames(study_config, study_config.subjects(study_config.current_subject).id);
fname = 'MergedDataBase_forBehavior.mat';
file2save = fullfile(N.searchFolder_2, fname);
if exist(file2save, 'file') && ~exist('MergedDataBase', 'var') && ~recompute
    load(file2save)
end

for subject_ind = subject_inds
    %subject_ind = 41;
    subject = study_config.subjects(subject_ind).id;
    study_config.current_subject = subject_ind;
    
    if strcmp(study_config.subjects(subject_ind).excluded, 'Yes')
        continue
    elseif exist('MergedDataBase', 'var') && ~recompute && sum(strcmp(MergedDataBase.ID, subject))>0
        continue
    else
        N = makeFolderFileNames(study_config, subject);
        filepath = N.searchFolder_2;
        
        % load DataBase
        load(fullfile(filepath, sprintf('%s_DataBase_withMocap.mat', subject)));
        
        if ~any(strcmp(fieldnames(Trials), 'ConfidenceLevel'))
            Trials = addvars(Trials, nan(size(Trials,1),1), 'NewVariableNames', {'ConfidenceLevel'});
        end
        
        if ~exist('MergedDataBase', 'var')
            MergedDataBase = removevars(Trials, {'urevent_seq', 'EEGComplete_trial', 'EEGComplete_details'});
        else
            MergedDataBase = [MergedDataBase; removevars(Trials, {'urevent_seq', 'EEGComplete_trial', 'EEGComplete_details'})];
        end
    end
end

%% save
if ~exist(file2save) || recompute
    save(file2save,'MergedDataBase')
end

%% Plot visibility for each trial
IDs = unique(MergedDataBase.ID);
trialsWithObject = strcmp(MergedDataBase.TrialType, 'WithObject');
trialsWithGoggles = strcmp(MergedDataBase.Condition, 'GogglesON');
for p = 1:numel(IDs)
    pLines = strcmp(MergedDataBase.ID, IDs{p});
    blocks = unique(MergedDataBase.Block(pLines & trialsWithGoggles));
    nBlocks = length(blocks);
    
    figure
    suptitle(IDs{p})
    for b = 1:nBlocks
        bLines = MergedDataBase.Block == blocks(b);
        trials = find(trialsWithObject & bLines & pLines);
        nTrials = length(trials);
        for t = 1:nTrials
            data = [MergedDataBase.PercTumblerVis_EC{trials(t)},MergedDataBase.PercTumblerVis_EO{trials(t)}];
            
            pl = subplot(nBlocks,5,t+(b-1)*5);
            hold on;
            plot(data, 'Marker', 'o', 'MarkerSize', 4, 'LineStyle', 'none');
                       xlim([0,1000]);
            ylim([0,100]);
            xline(length(MergedDataBase.PercTumblerVis_EC{trials(t)}), '--k', 'LineWidth', 2);
            yline(mean(MergedDataBase.PercTumblerVis_EO{trials(t)}), '--k', 'LineWidth', 1.25, 'Label', 'EO mean');
            
            switch MergedDataBase.Answer{trials(t)}
                case 'Absent'
                    c = 'w';
                case 'Present'
                    c = 'b';
            end
            
            plPos = get(pl, 'Position');
            diameter = 0.1;
            switch nBlocks
                case 3
                    ratio = 0.83;
                case 4
                    ratio = 0.65;
            end
            dims = [plPos(1)+(0.975-ratio*diameter)*plPos(3), plPos(2)+(0.975-diameter)*plPos(4),...
                ratio*diameter*plPos(3), diameter*plPos(4)];
            annotation('ellipse', dims, 'Color', 'b', 'FaceColor', c)
            
            if t == 1
                ylabel('Perc Visible');
            end
            
            if b == nBlocks
                xlabel('Sample Index');
            end
        end
    end
end

%% Plot visibility as histogram
visIndex = zeros(size(MergedDataBase,1),1);
perc_step = 5; edges = 0:perc_step:100;
tumbVis_EC = nan(size(MergedDataBase,1),length(edges)-1);
tumbVis_EO = nan(size(MergedDataBase,1),length(edges)-1);
for t = 1:size(MergedDataBase,1)
    if strcmp(MergedDataBase.TrialType{t}, 'WithObject')
        if strcmp(MergedDataBase.Condition{t}, 'GogglesON')
            visIndex(t) = sum(MergedDataBase.PercTumblerVis_EO{t}>0)/length(MergedDataBase.PercTumblerVis_EO{t});
            [tumbVis_EC(t,:),~] = histcounts(MergedDataBase.PercTumblerVis_EC{t}, edges,...
                'Normalization','probability');
            [tumbVis_EO(t,:),~] = histcounts(MergedDataBase.PercTumblerVis_EO{t}, edges,...
                'Normalization','probability');
        else
            visIndex(t) = 1;
        end
    end
end

plot_params.xlabel = true;
for subject_ind = subject_inds
    subject = study_config.subjects(subject_ind).id;
    figure;
    for p = 1:2
        subplot(1,2,p)
        if p == 1
            plot_params.title = 'Eyes Closed period';
            plot_params.ylabel = true;
            plot_params.color = [0 0.4470 0.7410];
            plotTumblerVisAcrossTrials(tumbVis_EC(strcmp(MergedDataBase.ID, subject),:), edges, plot_params);
        else
            plot_params.title = 'Eyes Open period';
            plot_params.ylabel = false;
            plot_params.color = [0.8500 0.3250 0.0980];
            plotTumblerVisAcrossTrials(tumbVis_EO(strcmp(MergedDataBase.ID, subject),:), edges, plot_params);
        end
    end
    suptitle(sprintf('%s - Mean proportion of tumbler visibility across trials (GogglesON, WithObject)',subject));
end


data = addvars(MergedDataBase,visIndex);
asw = zeros(size(data,1),1);
asw(strcmp(data.Answer,'Present')) = 1;
data.Answer = asw;
CorrectAnswer = zeros(size(data,1),1);
CorrectAnswer((data.Answer == 1 & strcmp(data.TrialType,'WithObject')) |...
    (data.Answer == 0 & strcmp(data.TrialType,'WithoutObject'))) = 1;
data = addvars(data,CorrectAnswer);
formula1 = ['Answer~1+Condition+TrialType+visIndex+',...
    'Condition:TrialType+(1+Condition+visIndex|ID)'];
% lme = fitlme(data, formula, 'DummyVarCoding', 'effects');
%
% [B,dev, stats] = mnrfit(visIndex,categorical(data.Answer))
% scatter(1:length(visIndex),sort(visIndex));

glme1 = fitglme(data, formula1, 'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'MPL',...
    'DummyVarCoding', 'effects', 'BinomialSize',1, 'DispersionFlag',false,...
    'CheckHessian', false, 'CovarianceMethod', 'conditional', 'CovariancePattern', 'FullCholesky');

anova(glme1)

formula2 = ['CorrectAnswer~1+Condition+TrialType+visIndex+',...
    'Condition:TrialType+(1+Condition+visIndex|ID)'];
glme2 = fitglme(data, formula2, 'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'MPL',...
    'DummyVarCoding', 'effects', 'BinomialSize',1, 'DispersionFlag',false,...
    'CheckHessian', false, 'CovarianceMethod', 'conditional', 'CovariancePattern', 'FullCholesky');

anova(glme2)

vect_ggON_Obj = zeros(1,glme2.NumCoefficients);vect_ggON_Obj(1) = 1;
vect_ggON_noObj = zeros(1,glme2.NumCoefficients);vect_ggON_noObj(1) = 1;
vect_ggOFF_Obj = zeros(1,glme2.NumCoefficients);vect_ggOFF_Obj(1) = 1;
vect_ggOFF_noObj = zeros(1,glme2.NumCoefficients);vect_ggOFF_noObj(1) = 1;
vect_ggON_Obj(2) = 1;vect_ggON_noObj(2) = 1;vect_ggOFF_Obj(2) = -1;vect_ggOFF_noObj(2) = -1;
vect_ggON_Obj(3) = 1;vect_ggON_noObj(3) = -1;vect_ggOFF_Obj(3) = 1;vect_ggOFF_noObj(3) = -1;
vect_ggON_Obj(5) = 1;vect_ggON_noObj(5) = -1;vect_ggOFF_Obj(5) = -1;vect_ggOFF_noObj(5) = 1;
[pval, Fval, DF1, DF2] = coefTest(glme2,vect_ggON_Obj-vect_ggOFF_Obj);
[pval, Fval, DF1, DF2] = coefTest(glme2,vect_ggON_noObj-vect_ggOFF_noObj);
[pval, Fval, DF1, DF2] = coefTest(glme2,vect_ggON_Obj-vect_ggON_noObj);
[pval, Fval, DF1, DF2] = coefTest(glme2,vect_ggOFF_Obj-vect_ggOFF_noObj);

alpha_crit = 0.05/4;
threshold = finv(1-alpha_crit, DF1, DF2);

estimates_link = sum([vect_ggOFF_noObj',vect_ggOFF_Obj',vect_ggON_noObj',vect_ggON_Obj'].*glme1.Coefficients.Estimate,1);
estimates = exp(estimates_link)./(1+exp(estimates_link));
figure;
scatter(estimates,[1,2,3,4]);
yticks([1,2,3,4]);
yticklabels({'GogglesOFF-NoObj','GogglesOFF-Obj','GogglesON-NoObj','GogglesON-Obj'})
xlabel('Answer estimate');

estimates_link = sum([vect_ggOFF_noObj',vect_ggOFF_Obj',vect_ggON_noObj',vect_ggON_Obj'].*glme2.Coefficients.Estimate,1);
estimates = exp(estimates_link)./(1+exp(estimates_link));
figure;
scatter(estimates,[1,2,3,4]);
yticks([1,2,3,4]);
yticklabels({'GogglesOFF-NoObj','GogglesOFF-Obj','GogglesON-NoObj','GogglesON-Obj'})
xlabel('Correct Answer estimate');


figure;
scatter(visIndex,asw+(rand(size(asw))*0.1-0.05))

