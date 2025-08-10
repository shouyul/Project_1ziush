function statsCV = computeCVStatistic(TrainingRes, TestingRes, RandomRes, features2test)
% From the results of the classification across CV partition, computes
% various statistical indicators to indicate performance for the tested
% parameters (here only the number of features to select)

uniqueLabels = unique(RandomRes.True_labels);
n_classes = length(uniqueLabels);

n_reps = size(RandomRes.Scores,2);
nb_folds = size(TrainingRes.Scores,2);
numberOfPointsInROC = 100;

len_feat=length(features2test);
Random_Err_allReps = nan(1,n_reps);
Random_Err_bal_allReps = nan(1,n_reps);
Train_Err = nan(len_feat, nb_folds);
Train_Err_bal = nan(len_feat, nb_folds);
Test_Err = nan(len_feat, nb_folds);
Test_Err_bal = nan(len_feat, nb_folds);

Random_CM = nan(n_reps, n_classes, n_classes);
Train_CM = nan(len_feat, nb_folds, n_classes, n_classes);
Test_CM = nan(len_feat, nb_folds, n_classes, n_classes);

if n_classes == 2
    %Random_ROC= cell(1,n_reps);
    Random_AUC = nan(1,n_reps);
    Train_ROC= cell(len_feat, nb_folds);
    Train_AUC = nan(len_feat, nb_folds);
    Test_ROC= cell(len_feat, nb_folds);
    Test_AUC = nan(len_feat, nb_folds);
else
    %Random_ROC= cell(n_classes,n_reps);
    Random_AUC = nan(n_classes,n_reps);
    %Train_ROC= cell(len_feat, nb_folds);
    Train_AUC = nan(n_classes, len_feat, nb_folds);
    %Test_ROC= cell(len_feat, nb_folds);
    Test_AUC = nan(n_classes, len_feat, nb_folds);
end

%% Random levels
for rep = 1:n_reps
    targets = RandomRes.True_labels;
    outputs = RandomRes.Scores(:,rep);
    if n_classes == 2
        outputs(outputs == -1) = 0;
        [Random_CM(rep,:,:), Random_Err_allReps(rep), Random_Err_bal_allReps(rep)] = computeConfusion([0,1], targets, outputs);
        
        %     [Random_Err(rep),Random_CM(rep,:,:),~,percentages] = confusion(...
        %         cat(1,~RandomRes.True_labels',RandomRes.True_labels'),cat(1,-RandomRes.Scores(:,rep)',RandomRes.Scores(:,rep)'));
        %     Random_Err_bal(rep)=mean(percentages(2,1:2));
        
        [tpr,fpr,~] = roc(cat(1,~RandomRes.True_labels',RandomRes.True_labels'),cat(1,-RandomRes.Scores(:,rep)',RandomRes.Scores(:,rep)'));
        % Uniformize ROC calculation
        %[Random_ROC{rep}.fpr, Random_ROC{rep}.tpr]=uniform_ROC(fpr{2}, tpr{2}, numberOfPointsInROC);
        %plot_roc(imposed_fpr, imposed_tpr)
        % AUC calculation
        Random_AUC(rep)=AUC_manual(fpr{2}, tpr{2});
    else
        [Random_CM(rep,:,:), Random_Err_allReps(rep), Random_Err_bal_allReps(rep)] = computeConfusion(uniqueLabels, targets, outputs);
        
        rocTargets = zeros(n_classes, length(RandomRes.True_labels));
        rocScores = zeros(n_classes, length(RandomRes.True_labels));
        for i = 1:length(RandomRes.True_labels)
            rocTargets(RandomRes.True_labels(i)+1,i) = 1;
            rocScores(RandomRes.Scores(i)+1,i) = 1;
        end
        [tpr,fpr,~] = roc(rocTargets,rocScores);
        
        for c = 1:n_classes
            Random_AUC(c,rep) = AUC_manual(fpr{c}, tpr{c});
        end
    end
end
Random_Err = mean(Random_Err_allReps);
Random_Err_signif = quantile(Random_Err_allReps, 0.05);
Random_Err_bal = mean(Random_Err_bal_allReps);
Random_Err_bal_signif = quantile(Random_Err_bal_allReps, 0.05);
%Random_CM = squeeze(mean(Random_CM,1));
Random_AUC = mean(Random_AUC,2);

for i=1:len_feat
    for fold=1:nb_folds
        if n_classes == 2
            %% Training
            targets = TrainingRes.True_labels{i,fold};
            outputs = TrainingRes.Scores{i,fold}(:,2)>0.5;
            [Train_CM(i,fold,:,:), Train_Err(i,fold), Train_Err_bal(i,fold)] = computeConfusion([0,1],targets, outputs);
            
            %         [Train_Err(i,fold),Train_CM(i,fold,:,:),~,percentages] = confusion(...
            %             cat(1,~TrainingRes.True_labels{i,fold}',TrainingRes.True_labels{i,fold}'),TrainingRes.Scores{i,fold}');
            %         Train_Err_bal(i,fold)=mean(percentages(2,1:2));
            
            [tpr,fpr,~] = roc(...
                cat(1,~TrainingRes.True_labels{i,fold}',TrainingRes.True_labels{i,fold}'),TrainingRes.Scores{i,fold}');
            % Uniformize ROC calculation
            [Train_ROC{i,fold}.fpr, Train_ROC{i,fold}.tpr]=uniform_ROC(fpr{2}, tpr{2}, numberOfPointsInROC);
            %plot_roc(imposed_fpr, imposed_tpr)
            % AUC calculation
            Train_AUC(i,fold)=AUC_manual(fpr{2}, tpr{2});
            
            %% Testing
            targets = TestingRes.True_labels{i,fold};
            outputs = TestingRes.Scores{i,fold}(:,2)>0.5;
            [Test_CM(i,fold,:,:), Test_Err(i,fold), Test_Err_bal(i,fold)] = computeConfusion([0,1],targets, outputs);
            
            %         [Test_Err(i,fold),Test_CM(i,fold,:,:),~,percentages] = confusion(...
            %             cat(1,~TestingRes.True_labels{i,fold}',TestingRes.True_labels{i,fold}'),TestingRes.Scores{i,fold}');
            %         Test_Err_bal(i,fold)=mean(percentages(2,1:2));
            
            [tpr,fpr,~] = roc(...
                cat(1,~TestingRes.True_labels{i,fold}',TestingRes.True_labels{i,fold}'),TestingRes.Scores{i,fold}');
            % Uniformize ROC calculation
            [Test_ROC{i,fold}.fpr, Test_ROC{i,fold}.tpr]=uniform_ROC(fpr{2}, tpr{2}, numberOfPointsInROC);
            %plot_roc(imposed_fpr, imposed_tpr)
            % AUC calculation
            Test_AUC(i,fold)=AUC_manual(fpr{2}, tpr{2});
            
            
            %{
        % classification error
        trainingErrors(i,fold) = classerror(TrainingRes.True_labels{i,fold}, TrainingRes.Scores{i,fold}(:,2)>=0.5);
        testErrors(i,fold) = classerror(TestingRes.True_labels{i,fold}, TestingRes.Scores{i,fold}(:,2)>=0.5);
        
        % confusion matrix
        trainingConfusion(i,fold,:,:)= confusionmat(TrainingRes.True_labels{i,fold}, TrainingRes.Scores{i,fold}(:,2)>=0.5);
        testConfusion(i,fold,:,:) = confusionmat(TestingRes.True_labels{i,fold}, TestingRes.Scores{i,fold}(:,2)>=0.5);
        
        % ROC and AUC
        trainingScores = TrainingRes.Scores{i,fold};
        trainingScoresForPosclass = trainingScores(:, 2);
        [X, Y, ~, AUC] = perfcurve(TrainingRes.True_labels{i,fold}, trainingScoresForPosclass, posclass,...
            'XVals', linspace(0,1,numberOfPointsInROC));
        %trainingROC{i,fold} = [X, Y];
        trainingAUC(i,fold)= AUC;
        
        testScores = TestingRes.Scores{i,fold};
        testScoresForPosclass = testScores(:,2);
        [X, Y, ~, AUC] = perfcurve(TestingRes.True_labels{i,fold}, testScoresForPosclass, posclass,...
            'XVals', linspace(0,1,numberOfPointsInROC));
        %testROC{i,fold} = [X, Y];
        testAUC(i,fold)= AUC;
            %}
        else
            %% Training
            targets = TrainingRes.True_labels{i,fold};
            outputs = TrainingRes.Scores{i,fold};
            [Train_CM(i,fold,:,:), Train_Err(i,fold), Train_Err_bal(i,fold)] = computeConfusion(uniqueLabels,targets, outputs);
            
            rocTargets = zeros(n_classes, length(targets));
            rocScores = zeros(n_classes, length(targets));
            for j = 1:length(targets)
                rocTargets(targets(j)+1,j) = 1;
                rocScores(outputs(j)+1,j) = 1;
            end
            [tpr,fpr,~] = roc(rocTargets,rocScores);
            
            % AUC calculation
            for c = 1:n_classes
                Train_AUC(c,i,fold) = AUC_manual(fpr{c}, tpr{c});
            end
            
            %% Testing
            targets = TestingRes.True_labels{i,fold};
            outputs = TestingRes.Scores{i,fold};
            [Test_CM(i,fold,:,:), Test_Err(i,fold), Test_Err_bal(i,fold)] = computeConfusion(uniqueLabels,targets, outputs);
            
            rocTargets = zeros(n_classes, length(targets));
            rocScores = zeros(n_classes, length(targets));
            for j = 1:length(targets)
                rocTargets(targets(j)+1,j) = 1;
                rocScores(outputs(j)+1,j) = 1;
            end
            [tpr,fpr,~] = roc(rocTargets,rocScores);
            
            % AUC calculation
            for c = 1:n_classes
                Test_AUC(c,i,fold) = AUC_manual(fpr{c}, tpr{c});
            end
        end
    end
end

%% output
statsCV.random.Err = Random_Err;
statsCV.random.Err_signif = Random_Err_signif;
statsCV.random.Err_bal = Random_Err_bal;
statsCV.random.Err_bal_signif = Random_Err_bal_signif;
%statsCV.random.CM = Random_CM;
statsCV.random.AUC = Random_AUC;

statsCV.training.Err = Train_Err;
statsCV.training.Err_bal = Train_Err_bal;
statsCV.training.CM = Train_CM;
if n_classes == 2
    statsCV.training.ROCs = Train_ROC;
end
statsCV.training.AUC = Train_AUC;

statsCV.testing.Err = Test_Err;
statsCV.testing.Err_bal = Test_Err_bal;
statsCV.testing.CM = Test_CM;
if n_classes == 2
    statsCV.testing.ROCs = Test_ROC;
end
statsCV.testing.AUC = Test_AUC;
end
