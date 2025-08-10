function DataBase2 = createTumbVisEvents(EEG, TumbVisChan, DataBase)
high_cutoff = 10;
thresh_vis_below = 10;
thresh_vis_above = 50;
bumper = EEG.srate/2;
wdw = 50;
thresh_vis_slope = 100*wdw/EEG.srate;

varNames = DataBase.Properties.VariableNames;
varTypes = varfun(@class, DataBase, 'OutputFormat', 'cell');
DataBase2 = table('Size', [0,8], 'VariableTypes', varTypes(1:8), 'VariableNames', varNames(1:8));
eventLats = [];
count = 0;

EEG_filt = custom_filter(EEG, [], high_cutoff);
for tr = 1:size(DataBase,1)
    if strcmp(DataBase.Condition{tr}, 'GogglesON') &&...
            strcmp(DataBase.TrialType{tr}, 'WithObject')
        latencies = ceil(EEG_filt.event(findInStructWithEmpties(EEG_filt.event, 'urevent', DataBase.urevent_seq(tr,2))).latency):...
            floor(EEG_filt.event(findInStructWithEmpties(EEG_filt.event, 'urevent', DataBase.urevent_seq(tr,3))).latency);
        
        figure;
        hold on
        p1 = plot(EEG_filt.data(TumbVisChan,latencies));
        yline(thresh_vis_below,'--k', 'Label', 'Vis. Threshold');
        above_th = EEG_filt.data(TumbVisChan,latencies) >= thresh_vis_below;
        %plot(find(above_th), EEG_filt.data(TumbVisChan,latencies(above_th)))
        differences = diff(EEG_filt.data(TumbVisChan,latencies));
        
        sum_diff_lats = 1+wdw/2:length(differences)-wdw/2;
        sum_differences = zeros(1,length(latencies));
        for l = sum_diff_lats
            sum_differences(l) = sum(differences(l-wdw/2:l+wdw/2));
        end
        %plot(differences+50)
        p2 = plot(sum_differences+50);
        yline(50+thresh_vis_slope,'--k', 'Label', 'Cum. slope Threshold');
        ylim([0,100])
        
        cand_events = find(diff(above_th)==1);
        if ~isempty(cand_events)
            for e = 1:length(cand_events)
                xl = xline(cand_events(e), '--r');
                leg = legend([p1,p2],{'Visbility',sprintf('Cumulative slope (%g samples window)',wdw)});
                if sum_differences(cand_events(e)) > thresh_vis_slope &&...
                        all(EEG_filt.data(TumbVisChan,latencies(cand_events(e))+(-bumper:0))<thresh_vis_above)&&...
                        any(EEG_filt.data(TumbVisChan,latencies(cand_events(e))+(1:bumper))>=thresh_vis_above)&&...
                        all(EEG_filt.data(TumbVisChan,latencies(cand_events(e))+(1:bumper))>=thresh_vis_below)
                    count = count+1;
                    DataBase2(count,:) = DataBase(tr,1:8);
                    eventLats = [eventLats; latencies(cand_events(e))+1];
%                     user_decision = input('Remove event? [0/1] ');
%                     if ~user_decision
%                         xline(cand_events(e), '--k');
%                     end
                else
%                     user_decision = input('Keep event? [0/1] ');                    
%                     if user_decision
%                         xline(cand_events(e), '--k');
%                     end
                end
                delete(xl)
                delete(leg)
            end
        end
        
        close gcf
    end
end

DataBase2 = addvars(DataBase2,eventLats,'NewVariableNames',{'TumbVisEventLat'});
end