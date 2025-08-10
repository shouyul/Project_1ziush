function plotERP(data, times_ms, mask, opes, opts)

%% options struct
if isempty(opts.ylims)
    ylims = [-inf,inf];
else
    ylims = opts.ylims;
end

%% Data operations
if opes.singleTrialNorm
    % single trial normalization
    singletrialbase = mean(data,2);
    ERPs = data./repmat(singletrialbase, [1,length(times_ms),1]);
else
    ERPs = data;
end

if opes.preStimBaseline
    preStimBase = mean(ERPs(:,times_ms<0,:),2);
    ERPs = ERPs./repmat(preStimBase, [1,size(ERPs,2),1]);
elseif ~isempty(opes.baseline)
    error('Check before using')
    % External baseline provided by the user
    siz = size(squeeze(opes.baseline));
    if siz == size(squeeze(ERPs(:,1,1)))
        ERPs = ERPs./repmat(opes.baseline, [1,size(ERPs,2),size(ERPs,3)]);
    else
        error('This type of baseline is no supported')
    end
end

% Average over all channels
ERPs = squeeze(mean(ERPs,1));

% Actual figure (no figure creation to allow subplotting)
hold on
switch opts.style
    case 'GrandAverage'
        % Average over trials
        trials = 1:size(ERPs,2);
        [meanERP, ERPstd] = computeERPbyCondition(ERPs, trials, opts.Smooth, opts.nSampsSmooth);
        line = plot(times_ms, meanERP', 'LineWidth',2);
        patch([times_ms, fliplr(times_ms)],[meanERP'-ERPstd', fliplr(meanERP'+ERPstd')],...
            line.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3);
        
%     case 'Phase-wise'
%         encoding_trials = contains(stc_time_rois.trialinfo.Phase,'Encoding');
%         
%         % Average over trials
%         [meanERP, ERPstd] = computeERPbyCondition(ERPs, encoding_trials, opts.Smooth, opts.nSampsSmooth);
%         line1 = plot(times_ms, meanERP, 'LineWidth',2);
%         patch([times_ms,fliplr(times_ms)],[meanERP-ERPstd, fliplr(meanERP+ERPstd)],...
%             line1.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
%         
%         % Average over trials
%         [meanERP, ERPstd] = computeERPbyCondition(ERPs, ~encoding_trials, opts.Smooth, opts.nSampsSmooth);
%         line2 = plot(times_ms, meanERP, 'LineWidth',2);
%         patch([times_ms,fliplr(times_ms)],[meanERP-ERPstd, fliplr(meanERP+ERPstd)],...
%             line2.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
%         
%         legend([line1,line2],{sprintf('Encoding (N=%d)',sum(encoding_trials)),...
%             sprintf('Test (N=%d)',sum(~encoding_trials))})
%         
%     case 'Condition-wise'
%         down_trials = contains(stc_time_rois.trialinfo.Condition,'DOWN');
%         up_trials = contains(stc_time_rois.trialinfo.Condition,'UP');
%         mix_trials = contains(stc_time_rois.trialinfo.Condition,'MIX');
%         
%         % Average over trials
%         [meanERP, ERPstd] = computeERPbyCondition(ERPs, down_trials, opts.Smooth, opts.nSampsSmooth);
%         line1 = plot(times_ms, meanERP, 'LineWidth',2);
%         patch([times_ms,fliplr(times_ms)],[meanERP-ERPstd, fliplr(meanERP+ERPstd)],...
%             line1.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
%         
%         % Average over trials
%         [meanERP, ERPstd] = computeERPbyCondition(ERPs, up_trials, opts.Smooth, opts.nSampsSmooth);
%         line2 = plot(times_ms, meanERP, 'LineWidth',2);
%         patch([times_ms,fliplr(times_ms)],[meanERP-ERPstd, fliplr(meanERP+ERPstd)],...
%             line2.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
%         
%         % Average over trials
%         [meanERP, ERPstd] = computeERPbyCondition(ERPs, mix_trials, opts.Smooth, opts.nSampsSmooth);
%         line3 = plot(times_ms, meanERP, 'LineWidth',2);
%         patch([times_ms,fliplr(times_ms)],[meanERP-ERPstd, fliplr(meanERP+ERPstd)],...
%             line3.Color, 'EdgeColor', 'none', 'FaceAlpha',0.3)
%         
%         legend([line1,line2,line3],{sprintf('Down (N=%d)',sum(down_trials)),...
%             sprintf('Up (N=%d)',sum(up_trials)), sprintf('Mix (N=%d)',sum(mix_trials))})
        
    otherwise
        error('Unknown Style')
end

yline(0,'--k');
xline(0,'--k');
%xline(0,'--k', 'Label', 'TumblerVisible');
ylim(ylims);

if ~isempty(mask)
    if any(mask)
        % make sure mask is a line
        if size(mask,1)>size(mask,2)
            mask = mask';
        end
        
        yl = ylim;
        yspan = diff(yl);
        intervals = mask2intervals(mask);
        for i = 1:size(intervals,1)
            inds = intervals(i,1):intervals(i,2);
            p = patch([times_ms(inds), fliplr(times_ms(inds))],...
                [repmat(yl(1)+yspan*0.05,1,length(inds)), repmat(yl(1)+yspan*0.1,1,length(inds))],...
                [0.5,0.5,0.5], 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
%         leg = [leg, p];
%         if numel(leg_desc) == 1
%             leg_desc = [leg_desc, 'Stat. signif. (against baseline)'];
%         else
%             leg_desc = [leg_desc, 'Stat. signif. (between plots)'];
%         end
    end
end

if opts.xlabel
    xlabel('Time (ms)');
end
if opts.ylabel
    ylabel('Potential (microV)');
end

    function [meanERP, ERPstd] = computeERPbyCondition(ERPdata, trials, smooth, n_samps_smooth)
        if smooth
            meanERP = squeeze(movmean(mean(ERPdata(:,trials),2), n_samps_smooth));
            ERPstd = squeeze(movmean(std(ERPdata(:,trials),[],2), n_samps_smooth));
        else
            meanERP = squeeze(mean(ERPdata(:,trials),2));
            ERPstd = squeeze(std(ERPdata(:,trials),[],2));
        end
    end
end