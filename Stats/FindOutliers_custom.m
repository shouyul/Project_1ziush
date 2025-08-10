function [outlier_ind] = FindOutliers_custom(input_vect, params, tail)
%FindOutliers_custom tries to automatically find outliers
%   Tests for normality by using SW test
%   If normal distributed uses modified z-score
%   If not normal distributed it uses adjusted boxplot
%
%   tail input allows to put the focus on the outlier whisker you want to
%   detect: 'both', 'left' (extreme low values), 'right' (extreme high values)
%
%   This function was updated on 25/11/2017 to take into account outliers
%   in the data when testing for normality. However, it will have further
%   updates in order to get a test that is robust enough
%
%   Modified from Janir Ramos da Cruz @ EPFL and IST

original_input = reshape(input_vect,1,length(input_vect)); % local copy
if any(isnan(original_input))
    input = original_input(~isnan(original_input));
else
    input = original_input;    
end

% Account for outliers in the normality testing
% Do modified z-score since it's robust, we don't care at this stage if
% skewed
ind_tmp = modified_zscore(input);
input_tmp = input; % local copy
% remove outliers for normality test
% but in other testing we use the original data
input_tmp(ind_tmp) = [];

[h,~,~] = swtest(input_tmp, 0.05);

if h==0
    disp('Distribution found to be normal, using zscore criterion')
    outlier_ind = modified_zscore(input);
elseif h==1
    disp('Distribution found to be skewed, using inner_fence criterion')
    outlier_ind = adjusted_Boxplot(input);
else
end

if any(isnan(original_input))
    original_inds = find(~isnan(original_input));
    outlier_ind = original_inds(outlier_ind);
end

% the function for the modified fold z-score
% ref:
% Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers",
% The ASQC Basic References in Quality Control: Statistical Techniques,
% Edward F. Mykytka, Ph.D., Editor.

    function [indices] = modified_zscore(X)
        % calculates the median absolute deviation
        X_mad = median(abs((X-repmat(median(X),1,length(X)))));
        % calculates the modified z-score
        X_zscore = 0.6745*(X-repmat(median(X),1,length(X)))./repmat(X_mad,1,length(X));
        % indices of the outliers
        switch tail
            case 'both'
                indices = find(abs(X_zscore) > params.z_criterion);
            case 'left'
                indices = find(X_zscore < -params.z_criterion);
            case 'right'
                indices = find(X_zscore > params.z_criterion);
            otherwise
                error('No such tail to find');
        end
    end

% the function for the adjusted boxplot
% ref:
% E. Vanderviere; M. Huber (2004). An Adjusted Boxplot for Skewed
% Distributions. COMPSTAT'2004 Symposium, Physica-Verlag/Springer.
%
% G. Brys; M. Hubert; P.J. Rousseeuw (2005). A Robustification of
% Independent Component Analysis. Journal of Chemometrics 19(5-7),
% pp. 364-375.
%
% J.W. Tukey (1977). Exploratory Data Analysis. Addison Wesley.

    function [indices] = adjusted_Boxplot(X)
        X_tmp = X; % local copy
        % computes the medcouple
        MC = medcouple(X_tmp);
        if abs(MC)>0.6
            disp('Highly skewed distribution, the detection of outliers may not be reliable')
        else
            disp(['Medcouple:', num2str(MC)])
        end
        
        skew_side_limit = params.skew_side_limit;
        opp_side_limit = params.opp_side_limit;       
        
        k1 = -opp_side_limit*(MC>=0) -skew_side_limit*(MC<0);
        k3 = skew_side_limit*(MC>=0) +opp_side_limit*(MC<0);
        % sort the data
        Y = sort(X_tmp);
        % compute the 25th percentile
        Q1 = median(Y(Y<median(Y)));
        % compute the 50th percentile
        Q2 = median(Y);
        % compute the 75th percentile
        Q3 = median(Y(Y>median(Y)));
        % compute the Interquartile Range
        IQR = Q3 - Q1;
        % find Q1 outliers
        Outlier_ind1 = find(Y < repmat(Q1 - params.inner_fence*exp(k1.*MC).*IQR,1,length(Y)));
        % find Q3 outliers
        Outlier_ind2 = find(Y > repmat(Q3 + params.inner_fence*exp(k3.*MC).*IQR,1,length(Y)));
        % The value of the outliers
        switch tail
            case 'both'
                Outliers = Y([Outlier_ind1 Outlier_ind2]);
            case 'left'
                Outliers = Y(Outlier_ind1);
            case 'right'
                Outliers = Y(Outlier_ind2);
            otherwise
                error('No such tail to find');
        end
        % Find the indice of the outliers
        Outliers_ind = ismember(X,Outliers);
        indices = find(Outliers_ind == 1);
    end
end