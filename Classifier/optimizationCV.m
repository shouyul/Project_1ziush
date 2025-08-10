function [Training, Testing, classifierInfo] = optimizationCV(features, labels, partition, orderMethod, features2test, model)
% Performs CV based on the given partition and the features presented.
% Fisher score ranking is used to gradually select the best features

uniqueLabels = unique(labels);
n_classes = length(uniqueLabels);

nb_folds = partition.NumTestSets;
len_feat=length(features2test);
trainingLabelTrue = cell(len_feat, nb_folds);
trainingPredictionScores = cell(len_feat, nb_folds);

testLabelTrue = cell(len_feat, nb_folds);
testPredictionScores = cell(len_feat, nb_folds);

%For Storing classifier infos
Nb_feat=reshape(repmat(features2test, [nb_folds,1]),[len_feat*nb_folds,1]);
Fold=repmat([1:nb_folds]', [len_feat, 1]);
%Sel_feats=cell(len_feat*nb_folds,1);
Sel_feats=cell(len_feat,nb_folds);
%Gamma=zeros(len_feat*nb_folds,1);
Gamma=zeros(len_feat,nb_folds);
%Delta=zeros(len_feat*nb_folds,1);
Delta=zeros(len_feat,nb_folds);

for i=1:len_feat
    feat=features2test(i);
    disp(['number of selected features:' , num2str(feat)])
    
    parfor fold=1:nb_folds
        % subfolds
        if isstruct(partition)
            trainingData = features(partition.training(:,fold), :);
            trainingLabels = labels(partition.training(:,fold));
            
            testData = features(partition.test(:,fold), :);
            testLabels = labels(partition.test(:,fold));
        else
            trainingData = features(partition.training(fold), :);
            trainingLabels = labels(partition.training(fold));
            
            testData = features(partition.test(fold), :);
            testLabels = labels(partition.test(fold));
        end
        % Fisher score for this training set
        [OrderInd, ~] = rankfeat(trainingData, trainingLabels, orderMethod);
        Selected_feat=OrderInd(1:feat);
        %Sel_feats{nb_folds*(i-1)+fold}=Selected_feat;
        Sel_feats{i,fold}=Selected_feat;
        
        trainingData_red = trainingData(:, Selected_feat);
        testData_red = testData(:, Selected_feat);
        
        % training
        switch model
            case {'svm', 'logistic'}
                if n_classes == 2
                    classifier = fitclinear(trainingData_red, trainingLabels, 'Learner', model);
                    % prediction
                    [~, trainingScores] = predict(classifier, trainingData_red);
                    [~, testScores] = predict(classifier, testData_red);
                    %Gamma(nb_folds*(i-1)+fold)=classifier.Gamma;
                    %Delta(nb_folds*(i-1)+fold)=classifier.Delta;
                else
                    classifier = fitcecoc(trainingData_red, trainingLabels, 'Learner', model);
                    % prediction
                    trainingScores = predict(classifier, trainingData_red);
                    testScores = predict(classifier, testData_red);
                    %Gamma(nb_folds*(i-1)+fold)=classifier.Gamma;
                    %Delta(nb_folds*(i-1)+fold)=classifier.Delta;
                end
                
            case {'linear', 'diaglinear'}
                %Options.MaxObjectiveEvaluations=20;
                %Options.ShowPlots=false;
                classifier = fitcdiscr(trainingData_red, trainingLabels, 'DiscrimType', model);
                % prediction
                [~, trainingScores] = predict(classifier, trainingData_red);
                [~, testScores] = predict(classifier, testData_red);
                %Gamma(nb_folds*(i-1)+fold)=classifier.Gamma;
                %Delta(nb_folds*(i-1)+fold)=classifier.Delta;
                Gamma(i,fold)=classifier.Gamma;
                Delta(i,fold)=classifier.Delta;
        end
        
        %{
        [train_tpr,train_fpr,train_thresholds] = roc(cat(1,~trainingLabels',trainingLabels'),trainingScores');
        [test_tpr,test_fpr,test_thresholds] = roc(cat(1,~testLabels',testLabels'),testScores');
        figure
        line([0,1],[0,1], 'LineStyle', '--', 'Color', 'k')
        hold on
        plot(train_fpr{2},train_tpr{2})
        plot(test_fpr{2},test_tpr{2})
        xlabel('False positive Rate')
        ylabel('True positive Rate')
        legend({'Chance level', 'ROC (Training)', 'ROC (Testing)'})
        title('Receiver operator characteristic')
        %}
        
        % store
        trainingLabelTrue{i,fold} = trainingLabels;
        trainingPredictionScores{i, fold} = trainingScores;
        testLabelTrue{i, fold} = testLabels;
        testPredictionScores{i, fold} = testScores;
    end
end

Training.True_labels=trainingLabelTrue;
Training.Scores=trainingPredictionScores;

Testing.True_labels=testLabelTrue;
Testing.Scores=testPredictionScores;

classifierInfo=table(Nb_feat, Fold,...
    reshape(Sel_feats',[len_feat*nb_folds,1]),...
    reshape(Gamma',[len_feat*nb_folds,1]),...
    reshape(Delta',[len_feat*nb_folds,1]));

% Training= struct('True_labels', trainingLabelTrue, ...
%     'Est_labels', trainingLabelEstimated, 'Scores', trainingPredictionScores);
% Testing= struct('True_labels', testLabelTrue, ...
%     'Est_labels', testLabelEstimated, 'Scores', testPredictionScores);
end

