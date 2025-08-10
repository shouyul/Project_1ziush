function [MC] = medcouple(X)
% Aims to compute the medcouple measure, a robust measure of skewness
% for a skewed distribution. Takes into account cases where the
% observations are equal to the median of the series.
%
% X should be a vector of observations (1 dimension)
%
% [MC] = medcouple(X) returns the following:
% MC    - vector with the medcouple measure of the data series
%
% Ref:
% G. Brys; M. Hubert; A. Struyf (2004). A Robust Measure of Skewness.
% Journal of Computational and Graphical Statistics 13(4), 996-1017.

X = reshape(X,length(X),1);
n = length(X);

s_X = sort(X);
X_med = median(X);
z = s_X - repmat(X_med,n,1);

ip = find(z(:,1)>=0); % These are the positions in z of z+
im = find(z(:,1)<=0); % These are the positions in z of z-
zp = z(ip,1); % z+ repeated to account for all cells in the matrix
zm = z(im,1); % z- repeated to account for all cells in the matrix

p = size(ip,1);
q = size(im,1);

[mi, mj] = ind2sub([p,q],1:p*q); % Positions of all combinations of z+
% and z- as elements in a pxq matrix

h = (zp(mi)+zm(mj))./(zp(mi)-zm(mj)); % same size as mi, mj

%% Special definition when z==0
ipz= find(zp==0); % row numbers of z+ = 0, i.e., X_{i} = median(X)
imz= find(zm==0); % row numbers of z- = 0, i.e., X_{i} = median(X)
piz = ismember(mi,ipz); % positions in mi where z+=0
pjz = ismember(mj,imz); % positions in mi where z-=0
zmember = piz+pjz; % same size as mi, mj
pijz = find(zmember == 2); % positions where z+ = z- = 0, i.e., X_{i} = X_{j} = median(X)
[indi,indj] = ind2sub([p,q],pijz); % pxq matrix position of the zero entries
indi = indi - min(indi) + 1; % row position of the zero entries as if they
% were in a separated matrix
indj = indj - min(indj) + 1; % column position of the zero entries as if they
% were in a separated matrix

for i=1:size(pijz,2)
    if (indi(i) + indj(i) - 1) > size(find(z==0),1)
        h(pijz(i)) = 1;
    elseif (indi(i) + indj(i) - 1) < size(find(z==0),1)
        h(pijz(i)) = -1;
    else
        h(pijz(i)) = 0;
    end
end

MC = median(h);
end

