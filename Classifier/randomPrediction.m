function RandomRes = randomPrediction(Labels_binary, n_reps)
% Random prediction
prior = sum(Labels_binary)/length(Labels_binary);
random_predict = rand(length(Labels_binary),n_reps);
random_predict(random_predict<=(1-prior)) = -1;
random_predict(random_predict>(1-prior)) = 1;
RandomRes.True_labels = Labels_binary;
RandomRes.Scores = random_predict;

Random_Err = nan(1,n_reps);
Random_Err_bal = nan(1,n_reps);
Random_CM = nan(n_reps, 2, 2);
Random_AUC = nan(1,n_reps);
for rep = 1:n_reps
    targets = Labels_binary;
    outputs = random_predict(:,rep);
    outputs(outputs == -1) = 0;
    [Random_CM(rep,:,:), Random_Err(rep), Random_Err_bal(rep)] = computeConfusion([0,1], targets, outputs);
    
    [tpr,fpr,~] = roc(cat(1,~RandomRes.True_labels',RandomRes.True_labels'),cat(1,-RandomRes.Scores(:,rep)',RandomRes.Scores(:,rep)'));
    % AUC calculation
    Random_AUC(rep)=AUC_manual(fpr{2}, tpr{2});
end
RandomRes.Err = Random_Err;
RandomRes.Err_bal = Random_Err_bal;
RandomRes.CM = Random_CM;
RandomRes.AUC = Random_AUC;
end