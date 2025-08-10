function CVstats = classifyWithFisher(Features, Labels, cfg, params)
if length(unique(Labels)) == 2
    % Binary classification
    %% Random classification:
    prior = sum(Labels)/length(Labels);
    random_predict = rand(length(Labels),1000);
    random_predict(random_predict<=(1-prior)) = -1;
    random_predict(random_predict>(1-prior)) = 1;
    RandomRes = struct();
    RandomRes.True_labels = Labels;
    RandomRes.Scores = random_predict;
else
    %% Random classification:
    uniqueLabels = unique(Labels);
    priors = nan(1,length(uniqueLabels));
    for u = 1:length(uniqueLabels)
        priors(u) = sum(Labels==uniqueLabels(u))/length(Labels);
    end
    
    random_predict = rand(length(Labels),1000);
    for u = 1:length(uniqueLabels)
        random_predict(random_predict<=sum(priors(1:u))) = uniqueLabels(u);
    end
    RandomRes = struct();
    RandomRes.True_labels = Labels;
    RandomRes.Scores = random_predict;
end

disp('Classification')
%ShrinkageAnalysis(Features_final, labels_target, OrderInd);
if strcmp(params.CV,'LOO')
    partition = cvpartition(Labels, 'LeaveOut');
elseif strcmp(params.CV,'LOOTR')
    Trials = unique(params.ChunksDependency);
    Labels_tr = zeros(size(Trials));
    for tr = 1:numel(Trials)
        if all(Labels(params.ChunksDependency == Trials(tr)))
            Labels_tr(tr) = 1;
        end
    end
    partition_tr = cvpartition(Labels_tr, 'LeaveOut');
    
    % create a custom struct for partition (not a cvpartition object) that
    % contains the necessary information for later
    partition.NumTestSets = numel(Trials);
    partition.training = true(numel(Labels),numel(Trials));
    partition.test = false(numel(Labels),numel(Trials));
    
    for fold = 1:numel(Trials)
        chunks_in_fold = params.ChunksDependency == Trials(partition_tr.test(fold));
        partition.training(chunks_in_fold,fold) = false;
        partition.test(chunks_in_fold,fold) = true;
    end
else
    partition = cvpartition(Labels, 'kfold', cfg.nbfolds);
end

% Mdl = fitcecoc(Features, Labels, 'Learners', cfg.model, 'CVPartition', partition, 'Options', statset('UseParallel',true));
% predictedLabels = kfoldPredict(Mdl,'Options', statset('UseParallel',true));
% figure
% confMat = confusionchart(Labels, predictedLabels,'RowSummary','row-normalized');
%
% [optMdl, HyperOptParams] = fitcecoc(Features, Labels, 'Learners', cfg.model,...
%     'OptimizeHyperparameters','auto', 'HyperparameterOptimizationOptions',...
%     struct('UseParallel', true, 'CVPartition', partition));
% predictedLabels_opt = resubPredict(optMdl,'Options', statset('UseParallel',true));
% figure
% confMat_opt = confusionchart(Labels, predictedLabels_opt,'RowSummary','row-normalized');

[TrainingRes, TestingRes, classifierInfo] = optimizationCV(Features, Labels, partition, 'fisher',...
    cfg.feat2test, cfg.model);
CVstats = computeCVStatistic(TrainingRes, TestingRes, RandomRes, cfg.feat2test);

if isfield(params,'phase')
    % Then it's Tumbler
    switch params.CV
        case 'folds'
            save(sprintf('%s%s_results-%s_%dfolds_%s', params.saveDataFolder, params.name, params.phase, cfg.nbfolds, params.suffix),...
                'TrainingRes', 'TestingRes', 'classifierInfo', 'CVstats');
        case 'LOO'
            save(sprintf('%s%s_results-%s_LOO_%s', params.saveDataFolder, params.name, params.phase, params.suffix),...
                'TrainingRes', 'TestingRes', 'classifierInfo', 'CVstats');
        case 'LOOTR'
            save(sprintf('%s%s_results-%s_LOOTR_%s', params.saveDataFolder, params.name, params.phase, params.suffix),...
                'TrainingRes', 'TestingRes', 'classifierInfo', 'CVstats');
    end
else
    % Then it's VEP
    switch params.CV
        case 'folds'
            save(sprintf('%s%s_results_%dfolds_%s', params.saveDataFolder, params.name, cfg.nbfolds, params.suffix),...
                'TrainingRes', 'TestingRes', 'classifierInfo', 'CVstats');
        case 'LOO'
            save(sprintf('%s%s_results_LOO_%s', params.saveDataFolder, params.name, params.suffix),...
                'TrainingRes', 'TestingRes', 'classifierInfo', 'CVstats');
        case 'LOOTR'
            save(sprintf('%s%s_results_LOOTR_%s', params.saveDataFolder, params.name, params.suffix),...
                'TrainingRes', 'TestingRes', 'classifierInfo', 'CVstats');
    end
end
end

