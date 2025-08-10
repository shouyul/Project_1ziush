function perc = percTumblerInCone(conePos, coneRM, tabRef, tumbWRTtab)
% figure;
% hold on
tic;

%% Compute tumbler shape & volume
tumbHei=0.06;
tumbDiam=0.06;

step_angle = 5; % in degrees
angles = step_angle:step_angle:360;
unitCircleX = nan(length(angles),1);
unitCircleY = nan(length(angles),1);
for a = 1:length(angles)
    unitCircleX(a) = cos(pi()*angles(a)/180);
    unitCircleY(a) = sin(pi()*angles(a)/180);
end

step_rad = tumbDiam/(2*5);
radii = step_rad:step_rad:tumbDiam/2;
bottomTumb = [0,0,0];
for r = radii
    bottomTumb = cat(1, bottomTumb,...
        [unitCircleX*r,unitCircleY*r,zeros(size(unitCircleX))]);
end
topTumb = bottomTumb;
topTumb(:,3) = tumbHei;

step_height = tumbHei/10;
heights = step_height:step_height:tumbHei-step_height;
% faceTumb = [];
sliceTumb = [];
for h = 1:length(heights)
    %     if mod(h,2) == 1
    %         select = 1:2:length(angles);
    %     else
    %         select = 2:2:length(angles);
    %     end
    %     faceTumb = cat(1, faceTumb,...
    %         [unitCircleX(select)*tumbDiam/2,unitCircleY(select)*tumbDiam/2,ones(size(unitCircleX(select)))*heights(h)]);
    newSlice = bottomTumb;
    newSlice(:,3) = heights(h);
    sliceTumb = cat(1, sliceTumb, newSlice);
end

%shapeTumb = cat(1,bottomTumb,faceTumb,topTumb);
volTumb = cat(1,bottomTumb,sliceTumb,topTumb);
% Translate to the right spot
for dim = 1:3
    %shapeTumb(:,d) = tabRef(d)+tumbWRTtab(d)+shapeTumb(:,d);
    volTumb(:,dim) = tabRef(dim)+tumbWRTtab(dim)+volTumb(:,dim);
end

% shp1 = alphaShape(shapeTumb(:,1),shapeTumb(:,2),shapeTumb(:,3));
% shp1.Alpha = 2.5;
% figure
% plot(shp1)


%% Compute cone shape
cone_length = 2; % in meters
cone_angle = 5*pi()/180; % in radians

cone = [0,0,0];
step_length = cone_length/20;
lengths = step_length:step_length:cone_length;
for l = 1:length(lengths)
    if mod(l,2) == 1
        select = 1:2:length(angles);
    else
        select = 2:2:length(angles);
    end
    newRad = lengths(l)*tan(cone_angle);
    cone = cat(1, cone,...
        [unitCircleX(select)*newRad,unitCircleY(select)*newRad,ones(size(unitCircleX(select)))*lengths(l)]);
end

n_pnts = size(coneRM,3);
if n_pnts == 1
    perc = percSinglePoint(cone, conePos, coneRM, volTumb);
else
    perc = nan(1,n_pnts);
    for i = 1:n_pnts
        perc(i) = percSinglePoint(cone, conePos(:,i), squeeze(coneRM(:,:,i)), volTumb);
    end
end

elapsedTime = seconds(toc);
elapsedTime.Format = 'mm:ss';
fprintf('Elapsed time for percTumblerInCone.m : %ss.\n', char(elapsedTime));

    function perc = percSinglePoint(cone, conePos, coneRM, volTumb)
        % Orient in the right direction
        for p = 1:size(cone,1)
            cone(p,:) = coneRM*cone(p,:)';
        end
        
        % Translate to the right spot
        for d = 1:3
            cone(:,d) = conePos(d)+cone(:,d);
        end
        
        coneShape = alphaShape(cone(:,1),cone(:,2),cone(:,3));
        coneShape.Alpha = 2.5;
        % plot(shp2)
        
        count = 0;
        for v = 1:size(volTumb,1)
            if inShape(coneShape, volTumb(v,1),volTumb(v,2),volTumb(v,3))
                count = count+1;
            end
        end
        perc = 100*count/size(volTumb,1);
    end
end