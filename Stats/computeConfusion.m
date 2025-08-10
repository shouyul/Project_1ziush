function [CM, err, err_bal] = computeConfusion(uniqueTargets, targets, outputs)

nTargets = length(uniqueTargets);

if nTargets == 2
    % target and output should be logical vectors
    % (1 for class of interest, 0 for the other)
    NPos_targets = sum(targets);
    NNeg_targets = length(targets)-NPos_targets;
    %     NPos_outputs = sum(outputs);
    %     NNeg_outputs = length(outputs)-NPos_outputs;
    Hits = sum(targets == 1 & outputs == 1);
    CRs = sum(targets == 0 & outputs == 0);
    Misses = sum(targets == 1 & outputs == 0);
    FAs = sum(targets == 0 & outputs == 1);
    err = (Misses + FAs)/length(targets);
    if NPos_targets > 0 && NNeg_targets > 0
        err_bal = 0.5*(Misses/NPos_targets + FAs/NNeg_targets);
    elseif NPos_targets == 0
        err_bal = FAs/NNeg_targets;
    else
        % NNeg_targets == 0
        err_bal = Misses/NPos_targets;
    end
    CM = [Hits,Misses;FAs,CRs];
else
    % Make sure targets and outputs are integer between 1 and nTargets
    truth = targets;
    prediction = outputs;
    for u = 1:nTargets
        truth(targets == uniqueTargets(u)) = u;
        prediction(outputs == uniqueTargets(u)) = u;
    end
    
    % Compute CM
    CM = zeros(nTargets, nTargets);
    for i = 1:length(truth)
        CM(truth(i),prediction(i)) = CM(truth(i),prediction(i))+1;
    end
    err = 1-(sum(diag(CM))/length(truth));
    err_bal = mean((sum(CM,2)-diag(CM))./sum(CM,2));
end
end