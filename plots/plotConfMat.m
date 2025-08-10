function plotConfMat(varargin)
%PLOTCONFMAT plots the confusion matrix with colorscale, absolute numbers
%   and precision normalized percentages
%
%   usage:
%   PLOTCONFMAT(confmat) plots the confmat with integers 1 to n as class labels
%   PLOTCONFMAT(confmat, labels) plots the confmat with the specified labels
%
%   Vahe Tshitoyan
%   20/08/2017
%
%   Arguments
%   confmat:                a square confusion matrix
%   labels (optional):      vector of class labels
%   percOnly (optional):    display the repartition percentage only.
%                           If set to false, the number of trials
%                           falling into each category will be displayed too.
%                           default: false

% number of arguments
switch (nargin)
    case 0
        confmat = 1;
        labels = {'1'};
        percOnly = false;
    case 1
        confmat = varargin{1};
        labels = 1:size(confmat, 1);
        percOnly = false;
    case 2
        confmat = varargin{1};
        labels = varargin{2};
        percOnly = false;
    otherwise
        confmat = varargin{1};
        labels = varargin{2};
        percOnly = varargin{3};
        
end

confmat(isnan(confmat))=0; % in case there are NaN elements
numlabels = size(confmat, 1); % number of labels

% calculate the percentage accuracies
confpercent = 100*confmat./repmat(sum(confmat, 1),numlabels,1);

% plotting the colors
imagesc(confpercent);
title(sprintf('Accuracy: %.2f%%', 100*trace(confmat)/sum(confmat(:))));
ylabel('Answer'); xlabel('Truth');

% set the colormap
colormap(flipud(gray));

% Create strings from the matrix values and remove spaces
if percOnly
    textStrings = num2str(confpercent(:), '%.1f%%\n');
else
    textStrings = num2str([confpercent(:), confmat(:)], '%.1f%%\n%.1f\n');
end
textStrings = strtrim(cellstr(textStrings));

% Create x and y coordinates for the strings and plot them
[x,y] = meshgrid(1:numlabels);
hStrings = text(x(:),y(:),textStrings(:), ...
    'HorizontalAlignment','center');

% Get the middle value of the color range
midValue = mean(get(gca,'CLim'));

% Choose white or black for the text color of the strings so
% they can be easily seen over the background color
textColors = repmat(confpercent(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2));

% Setting the axis labels
set(gca,'XTick',1:numlabels,...
    'XTickLabel',labels,...
    'YTick',1:numlabels,...
    'YTickLabel',labels,...
    'YTickLabelRotation',90,...
    'TickLength',[0 0]);
end