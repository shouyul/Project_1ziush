function [EEG] = IC_categorization(EEG, thresholds, user, type2inspect)
% Attributes classes from a compositional label array outputed by ICLabel
% Use the thresholds in input to determine possible categories for each
% component. Any case is possible: ICs ending up with 0, 1 or more
% attributed classes.
%
% Inputs:
%   EEG                     - EEG struct from eeglab where iclabel has been performed
%   thresholds              - thresholds above which a class is attributed to the IC [1 x nb_classes]
%   user                    - (optional) [true false] Ask the user if the labelling is undecided.
%                               default is false.
%   type2inspect            - (optional) ['all' 'brain' 'muscle' 'eye' 'heart' 'line' 'chan' 'other']
%                               In case 'user' is set to true:
%                               Which components should be inspected depending on the class labels that are eligible
%                               ex: 'brain' disentangle only components where brain is eligible
%                               default is 'brain'
%
% Outputs (written in the EEG.etc struct):
%   mostProbableClass       - vector of the most probable class per IC
%                              (the prediction value still has to be above the threshold)
%   eligible_classes        - vector of the number of eligible classes per IC
%   reportTable             - char array reporting the categorization of
%                               ICs (which classes were above the threshold and in which order)

classes = EEG.etc.ic_classification.ICLabel.classes;
classifications = EEG.etc.ic_classification.ICLabel.classifications;

mostProbableClass = zeros(size(classifications,1),1);
eligible_classes = [];
n_eligible_classes = zeros(size(classifications,1),1);
reportTable=[];

for IC=1:size(classifications,1)
    
    % Find predictions above threshold
    above_th = classifications(IC,:) >= thresholds;
    if any(above_th)
    n_eligible_classes(IC) = sum(above_th);
    
    % Sort the classes by prediction percentage    
    [vals, order] = sort(classifications(IC, above_th), 'descend');
    inds_unsorted = find(above_th);
    inds = inds_unsorted(order);
    
    % Fill report variables
    reportString = '';
    for c = 1:length(vals)
        eligible_classes(IC, c) = inds(c);
        if c == 1
            mostProbableClass(IC) = inds(c);
        end
        
        if strcmp(reportString, '')
            reportString = [num2str(c), '. ', classes{inds(c)},...
                ' (', num2str(vals(c)*100, '%.2f'),'% >= ', num2str(thresholds(inds(c))*100, '%d'), '%)'];
        else
            reportString = [reportString, ' ', num2str(c), '. ', classes{inds(c)},...
                ' (', num2str(vals(c)*100,'%.2f'),'% >= ', num2str(thresholds(inds(c))*100,'%d'), '%)'];
        end
    end
    
    
%     i = 1;
%     val=vals(i); ind= inds(i);
%     while (i<= length(classes) && val >= thresholds(ind))
%         eligible_classes(IC, i) = ind;
%         
%         if i==1
%             mostProbableClass(IC) = ind;
%         end
%         
%         if strcmp(reportString, '')
%             reportString = [num2str(i), '. ', classes{ind},...
%                 ' (', num2str(val*100, '%.2f'),'% >= ', num2str(thresholds(ind)*100, '%d'), '%)'];
%         else
%             reportString = [reportString, ' ', num2str(i), '. ', classes{ind},...
%                 ' (', num2str(val*100,'%.2f'),'% >= ', num2str(thresholds(ind)*100,'%d'), '%)'];
%         end
%         
%         i = i+1;
%         val=vals(i); ind= inds(i);
%     end
    
    reportString = string(reportString);
    reportTable = char([reportTable; reportString]);
    else
        warning('No class above threshold for IC%d', IC)
        eligible_classes(IC, 1) = 0; % Important to keep eligible classes with the right length        
    end
end

%% User defined labels (when undecided) - optional
if ~exist('user', 'var')
    user = false;
    type2inspect = 'none';
elseif ~exist('type2inspect', 'var')
    type2inspect = 'brain';
end

if user
    userSelectedClass = mostProbableClass;
    tmp_el_classes = eligible_classes;
    tmp_n_el_classes = n_eligible_classes;
    
    % automatic rejection of brain ICs when their associated dipole has too high rv
    % and when they are located outside the brain
    if isfield(EEG, 'dipfit') && ~isempty(EEG.dipfit)
        brain_candidates = find(sum(eligible_classes == 1,2));
        
        for comp = brain_candidates'
            if EEG.dipfit.model(comp).rv > 0.15 || sqrt(sum(EEG.dipfit.model(comp).posxyz.^2,2)) > 85 %85 is the standard head radius
                if n_eligible_classes(comp) == 1
                    % no other class passed the threshold
                    userSelectedClass(comp) = 7; % classify as other
                else% other classes passed the threshold
                    brain_ind = find(eligible_classes(comp,:)==1);
                    if brain_ind==1
                        % Brain was the most probable
                        userSelectedClass(comp) = eligible_classes(comp,2);
                        tmp_el_classes(comp, brain_ind:end-1) = eligible_classes(comp,brain_ind+1:end);
                    elseif brain_ind==length(eligible_classes(comp,:))
                        tmp_el_classes(comp, end) = 0;
                    else
                        tmp_el_classes(comp, brain_ind:end-1) = eligible_classes(comp,brain_ind+1:end);
                    end
                    tmp_n_el_classes(comp) = tmp_n_el_classes(comp)-1;
                end
            end
        end
    end
    
    switch lower(type2inspect)
        case 'all'
            undecided = tmp_n_el_classes>1;
        case 'brain'
            undecided = tmp_n_el_classes>1 & any(tmp_el_classes == 1,2);
        case 'muscle'
            undecided = tmp_n_el_classes>1 & any(tmp_el_classes == 2,2);
        case 'eye'
            undecided = tmp_n_el_classes>1 & any(tmp_el_classes == 3,2);
        case 'heart'
            undecided = tmp_n_el_classes>1 & any(tmp_el_classes == 4,2);
        case 'line'
            undecided = tmp_n_el_classes>1 & any(tmp_el_classes == 5,2);
        case 'chan'
            undecided = tmp_n_el_classes>1 & any(tmp_el_classes == 6,2);
        case 'other'
            undecided = tmp_n_el_classes>1 & any(tmp_el_classes == 7,2);
        otherwise
            error('Not a valid class')
    end
    
    % Add comps where 'Other' is the only class
    %undecided = undecided | (tmp_n_el_classes == 1 & any(tmp_el_classes == 7,2));
    
    % Add comps where no classes pass threshold
    undecided = undecided | tmp_n_el_classes == 0;
    
    % Convert to indices
    undecided = find(undecided)';
    
    freq_range = [0.1 60]; % in Hz
    spec_opt = {'freqrange', freq_range}; % cell array of options which are passed to spectopo()
    erp_opt = {}; % cell array of options which are passed to erpimage()
    for comp = undecided
        %rmpath('C:\Program Files\MATLAB\R2019a\toolbox\matlab\uicomponents\uicomponents\+matlab\+ui\+internal\+mixin\@Scrollable');
        pop_prop_extended(EEG, 0, comp, NaN, spec_opt , erp_opt, 1, 'ICLabel');
        disp('Brain:1 Muscle:2 Eye:3 Heart:4 LN:5 ChanNoise:6 Other:7');
        userSelectedClass(comp) = input('Which label to assign ? ');
        close(gcf)
        %addpath('C:\Program Files\MATLAB\R2019a\toolbox\matlab\uicomponents\uicomponents\+matlab\+ui\+internal\+mixin\@Scrollable','-end');
    end
else
    userSelectedClass = zeros(size(classifications,1),1);
end

%% Store everything
EEG.etc.ic_classification.ICLabel.detectionThresholds = thresholds;
EEG.etc.ic_classification.ICLabel.mostProbableClass = mostProbableClass;
EEG.etc.ic_classification.ICLabel.eligibleClasses = eligible_classes; % ordered from the most to the least probable
EEG.etc.ic_classification.ICLabel.sizeEligibleClasses = n_eligible_classes;
EEG.etc.ic_classification.ICLabel.report = reportTable;
EEG.etc.ic_classification.ICLabel.userSelectedClass = userSelectedClass;
EEG.etc.ic_classification.ICLabel.userInspectedType = lower(type2inspect);
end