function [bi_means, bi_STDs] = Biweight_custom(data, N_samples, censor_param)
% Biweight_custom calculates the biweight mean and std
%   data - is the input matrix, where the rows are the variables and
%       columns are the samples
%   outputs:
%   bi_mean - the biweight mean
%   bi_std - the biweight std
%
%   Ref:
%   Hoaglin, D., Mosteller, F., and J. Tukey 1983. Understanding Robust
%       and Exploratory Data Analysis, John Wiley and Sons, New York, 447 pp.
%   Lanzante JR. 1996. Resistant, robust and non-parametric techniques for
%       the analysis of climate data: theory and examples, including applications
%       to historical radiosonde station data. International Journal of
%       Climatology, vol.16, 1197-1226.

% Script adapted from Da Cruz

Medians = median(data, 2); % median for each channel
MADs = median(abs(data-Medians),2); % median absolute deviation that is the median of the sample
% of the absolute values of the differences from the median
weights = (data-Medians)./(censor_param*MADs); % weights for the computation of the biweight mean
weights(abs(weights)>1) = 1; % censoring of the weights

% Useful quantities (calculated on weights matrices):
Aw = (1-weights.^2).^2;
Bw = (1-weights.^2).^4;
Cw = (1-weights.^2);
Dw = (1-5*weights.^2);

bi_means = Medians + sum((data-Medians).*Aw,2)./sum(Aw,2); % computation of biwheight mean
bi_STDs = sqrt(N_samples*sum((data-Medians).^2.*Bw,2))./abs(sum(Cw.*Dw,2)); % computation of biwheight std
end