function doFisherAnalysis(Features, Labels, chan_labels, dim2_labels, max_feat, params)
disp('Fisher Analysis')
[OrderInd, scoreOfFeatures] = rankfeat(Features, Labels, 'fisher');

if isfield(params, 'phase')
    % It's Tumbler
    visualize_Fischer(OrderInd, scoreOfFeatures, 'all', chan_labels, dim2_labels);
    saveCurrentFig(params.saveFigFolder, sprintf('%s_FScore_allFeats-%s_%s', params.name, params.phase, params.suffix),...
        {'png','svg'}, [1000,1000])
    visualize_Fischer(OrderInd, scoreOfFeatures, max_feat, chan_labels, dim2_labels);
    saveCurrentFig(params.saveFigFolder, sprintf('%s_FScore_%dbestFeats-%s_%s', params.name, max_feat, params.phase, params.suffix),...
        {'png','svg'}, [1000,1000])
else
    % It's VEP
    visualize_Fischer(OrderInd, scoreOfFeatures, 'all', chan_labels, dim2_labels);
    saveCurrentFig(params.saveFigFolder, sprintf('%s_FScore_allFeats_%s', params.name, params.suffix),...
        {'png'}, [1000,1000])
    visualize_Fischer(OrderInd, scoreOfFeatures, max_feat, chan_labels, dim2_labels);
    saveCurrentFig(params.saveFigFolder, sprintf('%s_FScore_%dbestFeats_%s', params.name, max_feat, params.suffix),...
        {'png'}, [1000,1000])
end
end