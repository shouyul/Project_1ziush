function plotTumblerVisAcrossTrials(tumbVis, edges, params)
bar(edges(1:end-1)+diff(edges)/2,mean(tumbVis,1, 'omitnan'),1,...
    'EdgeColor', params.color, 'FaceColor', params.color);
ylim([0,1]);
if params.xlabel
    xlabel('% of tumbler visible');
end
if params.ylabel
    ylabel('Mean proportion in trials');
end
title(params.title);
end