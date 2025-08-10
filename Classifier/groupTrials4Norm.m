function tags = groupTrials4Norm(trialInfo, trialsGroup)
tags = zeros(size(trialInfo,1),1);
for s = 1:numel(trialsGroup)
    switch trialsGroup{s}
        case 'perSubject'
            subjs = unique(trialInfo.ID);
            for a = 1:numel(subjs)
                tags(strcmp(trialInfo.ID,subjs{a})) = tags(strcmp(trialInfo.ID,subjs{a}))+a;
            end
        case 'perCondition'
            conds = unique(trialInfo.Condition);
            for a = 1:numel(conds)
                tags(strcmp(trialInfo.Condition,conds{a})) = tags(strcmp(trialInfo.Condition,conds{a}))+10*a;
            end
        case 'perBlock'
            blks = unique(trialInfo.Block);
            for a = 1:length(blks)
                tags(trialInfo.Block==blks(a)) = tags(trialInfo.Block==blks(a))+100*a;
            end
        case 'perTrialType'
            TTs = unique(trialInfo.TrialType);
            for a = 1:numel(TTs)
                tags(strcmp(trialInfo.TrialType,TTs{a})) = tags(strcmp(trialInfo.TrialType,TTs{a}))+1000*a;
            end
        case 'perAnswer'
            Aswrs = unique(trialInfo.Answer);
            for a = 1:numel(Aswrs)
                tags(strcmp(trialInfo.Answer,Aswrs{a})) = tags(strcmp(trialInfo.Answer,Aswrs{a}))+10000*a;
            end
        otherwise
            error('Unknown trials grouping')
    end
end
end