function visualize_Fischer(OrderInd,scoreOfFeatures,max_feat,chan_labels,dim2_labels)
% Plot the Fischer score in 2D image to visualize which features are the most discriminant
% Inputs:
% OrderInd          - Vector of features' ranking as outputted by rankfeat()
% scoreOfFeatures   - Vector of features' score (corresponding to
%                       OrderInd). orderedPower ouput from rankfeat().
% max_feat          - Maximum number of features to plot (if not 'all', the
%                       features excluded will automatically appear with a 0 score on the plot.
% chan_labels       - labels of channels ({EEG.chanlocs.labels})
% dim2_labels       - labels of second dimension (to write on the ticks)

nb_dim2 = numel(dim2_labels);
nb_chans = numel(chan_labels);
[Chan_inds,dim2_inds] = ind2sub([nb_chans, nb_dim2], OrderInd);

switch max_feat
    case 'all'
        max_feat=numel(scoreOfFeatures);
        TITLE=['Visualisation of all features according to fisher score'];
    otherwise
        TITLE=['Visualisation of the ', num2str(max_feat),' best features according to fisher score'];
end

data2plot=zeros(nb_chans, nb_dim2);
for i=1:max_feat
    data2plot(Chan_inds(i), dim2_inds(i))=scoreOfFeatures(i);
end

if nb_chans > 10
    for i=1:nb_chans
        if rem(i,2)==0
            chan_labels{i}=[chan_labels{i},'-------'];
            %     elseif rem(i,3)==0
            %         chan_labels{i}=[chan_labels{i},'-------','-------'];
        end
    end
end

figure
imagesc(data2plot, [0,0.5])
% Hardcoded
if nb_dim2 <= 100
    xlabel('Frequency [Hz]')
else
    xlabel('Time [Samples]')
end
if nb_dim2 <= 5
    xticks(1:nb_dim2)
    xticklabels(dim2_labels)
elseif nb_dim2 <= 20
    xticks(1:2:nb_dim2)
    xticklabels(dim2_labels(1:2:end))
elseif nb_dim2 <= 100
    xticks(1:5:nb_dim2)
    xticklabels(dim2_labels(1:5:end))
else
    xticks(1:100:nb_dim2)
    xticklabels(dim2_labels(1:100:end))
    xtickangle(45)
end
ylabel('Channel')
yticks(1:nb_chans)
yticklabels(chan_labels)
c=colorbar();
c.Label.String=sprintf('Fisher score, max = %.2f', max(data2plot, [], 'all'));
c.Label.FontSize=14;
title(TITLE)
set(gca, 'Fontsize', 12)
end