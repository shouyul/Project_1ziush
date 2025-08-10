function intervals = getIntervals(EEG, intervalType, timeMargin, mergeOverlaps)

events = EEG.event;
Trials = EEG.etc.TrialsInspection;
validTrials = find([Trials.Trial_perc] == 100);

intervals = [];
inter_ind = 1;

for tr = 1:length(validTrials)
    blk = findInStructWithEmpties(events, 'block', Trials.BlockInd(validTrials(tr)));
    trl = findInStructWithEmpties(events, 'trial', Trials.TrialInd(validTrials(tr)));
    
    trial_evts = intersect(blk,trl);
    
    switch lower(intervalType)
        case 'fulltrial'
            % Keep only trials for which the time margin value doesn't include Nans
            if timeMargin > 0
                if Trials.MaxBufferPre(validTrials(tr)) <= timeMargin+(1/EEG.srate) || Trials.MaxBufferPost(validTrials(tr)) <= timeMargin+(1/EEG.srate)
                    % Add (1/EEG.srate) margin to avoid bad surprises                   
                    warning('Buffer too long: Had to remove trial %d in block %d', Trials.TrialInd(validTrials(tr)), Trials.BlockInd(validTrials(tr)));
                    continue
                end
            end
            
            trial_start_ev = [events.urevent] == Trials.Trial_urevent_seq(validTrials(tr),1);            
            intervals(inter_ind,1) = EEG.times(round(events(trial_start_ev).latency))/1000;
            
            trial_end_ev = [events.urevent] == Trials.Trial_urevent_seq(validTrials(tr),end);            
            intervals(inter_ind,2) = EEG.times(round(events(trial_end_ev).latency))/1000;           
            
%             trial_start_ev = intersect(find(strcmp({events(:).type},'TrialStart')),trial_evts);
%             if isempty(trial_start_ev)
%                 disp('problem')
%             else
%                 intervals(inter_ind,1) = EEG.times(round(events(trial_start_ev).latency))/1000;
%             end
%             
%             trial_end_ev = intersect(find(strcmp({events(:).type},'TrialEnd')),trial_evts);
%             if isempty(trial_end_ev)
%                 disp('problem')
%             else
%                 intervals(inter_ind,2) = EEG.times(round(events(trial_end_ev).latency))/1000;
%             end
            inter_ind = inter_ind+1;          
    end
end

if timeMargin > 0
    intervals(:,1) = intervals(:,1)-timeMargin;
    intervals(:,2) = intervals(:,2)+timeMargin;
end

if mergeOverlaps
    j=1;
    while j < size(intervals,1)
        if intervals(j,2)>=intervals(j+1,1)
            % Merge intervals
            intervals(j,2) = intervals(j+1,2);
            % Delete i+1 interval
            intervals(j+1,:) = [];
        else
            j=j+1;
        end
    end
end
end