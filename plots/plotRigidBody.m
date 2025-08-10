function plotRigidBody(RB, latencies, plot_style, plot_cone, tabRef, tumbWRTtab)
if ischar(plot_style)
    switch plot_style
        case 'all'
            % Plot all latencies as individual points
            step_samp = 1;
        case 'ave'
            % Plot the average over all latencies
            
    end
else
    % Number interpreted as a step
    step_samp = plot_style;
end

hold on
%% Plot table
tabLen=1.2;
tabWid=0.8;
tabDep=0.02;
[Xtab, Ytab] = meshgrid(tabRef(1):0.1:tabRef(1)+tabLen,tabRef(2):0.1:tabRef(2)+tabWid);
Ztab = ones(size(Xtab))*tabRef(3);
surf(Xtab,Ytab,Ztab, 'LineStyle', 'none', 'FaceAlpha', 1, 'FaceColor', [0.94,0.87,0.80]);

%% Plot tumbler
tumbHei=0.06;
tumbDiam=0.06;
[Xtumb,Ytumb,Ztumb] = cylinder(tumbDiam/2,100); % cylinder of length 1m
trans1=makehgtform('translate',[tabRef(1)+tumbWRTtab(1),tabRef(2)+tumbWRTtab(2),tabRef(3)+tumbWRTtab(3)],...
    'scale', [1,1,tumbHei]);
surf(Xtumb,Ytumb,Ztumb, 'Parent', hgtransform('Matrix',trans1),...
    'LineStyle', 'none', 'FaceAlpha', 1, 'FaceColor', [0.5,0.5,0.5]);
% cap for the tumbler
patch('XData',tabRef(1)+tumbWRTtab(1)+Xtumb(2,:),...
    'YData',tabRef(2)+tumbWRTtab(2)+Ytumb(2,:),...
    'ZData',tabRef(3)+tumbWRTtab(3)+Ztumb(2,:)*tumbHei,...
    'LineStyle', 'none', 'FaceAlpha', 1, 'FaceColor', [0.5,0.5,0.5]);

%% Plot RB and cone
if plot_cone
    % Compute Rotation matrices for cone plotting
    cone_RM = cat(2,reshape(RB.Zvector,size(RB.pos,1),1,size(RB.pos,2)),...
        reshape(RB.Xvector,size(RB.pos,1),1,size(RB.pos,2)),...
        reshape(RB.Yvector,size(RB.pos,1),1,size(RB.pos,2)));
end

% Scale vector length (default = 1m)
RB.Xvector = RB.Xvector.*0.1;
RB.Yvector = RB.Yvector.*0.25;
RB.Zvector = RB.Zvector.*0.1;

if exist('step_samp','var')
    for lat = latencies(1:step_samp:end)
        plotSinglePoint(RB.pos(:,lat), RB.Xvector(:,lat), RB.Yvector(:,lat), RB.Zvector(:,lat), squeeze(cone_RM(:,:,lat)));
        %     vector_x = drawArrow3([RB.pos(1,lat),RB.pos(1,lat)+RB.Xvector(1,lat)],...
        %         [RB.pos(2,lat),RB.pos(2,lat)+RB.Xvector(2,lat)],...
        %         [RB.pos(3,lat),RB.pos(3,lat)+RB.Xvector(3,lat)],...
        %         'linewidth',2.5,'color','b');
        %     vector_y = drawArrow3([RB.pos(1,lat),RB.pos(1,lat)+RB.Yvector(1,lat)],...
        %         [RB.pos(2,lat),RB.pos(2,lat)+RB.Yvector(2,lat)],...
        %         [RB.pos(3,lat),RB.pos(3,lat)+RB.Yvector(3,lat)],...
        %         'linewidth',2.5,'color','r');
        %     vector_z = drawArrow3([RB.pos(1,lat),RB.pos(1,lat)+RB.Zvector(1,lat)],...
        %         [RB.pos(2,lat),RB.pos(2,lat)+RB.Zvector(2,lat)],...
        %         [RB.pos(3,lat),RB.pos(3,lat)+RB.Zvector(3,lat)],...
        %         'linewidth',2.5,'color','g');
        %
        %     if plot_cone
        %         cone_length = 1;
        %         cone_angle = 5*pi()/180; % in radians
        %         max_radius = tan(cone_angle);
        %         [ax, ang] = rotMat2axisAngle(squeeze(cone_RM(:,:,lat)));
        %
        %         [Xcone,Ycone,Zcone] = cylinder([0,max_radius],100); % cone of length 1
        %         trans=makehgtform('scale', [1,1,cone_length],...
        %             'translate',[RB.pos(1,lat),RB.pos(2,lat),RB.pos(3,lat)],...
        %             'axisrotate', ax, ang);
        %         surf(Xcone,Ycone,Zcone, 'Parent', hgtransform('Matrix',trans),...
        %             'LineStyle', 'none', 'FaceAlpha', 0.25);
        %     end
    end
else
    plotSinglePoint(mean(RB.pos(:,latencies),2),...
        mean(RB.Xvector(:,latencies),2),...
        mean(RB.Yvector(:,latencies),2),...
        mean(RB.Zvector(:,latencies),2),...
        squeeze(mean(cone_RM(:,:,latencies),3)));
end
grid on
xlabel('-X');
ylabel('Z');
zlabel('Y');
axis equal
xlim([tabRef(1),tabRef(1)+tabLen]);
ylim([-inf,tabRef(2)+tabWid]);
zlim([tabRef(3),inf]);
view(-125,15);

    function plotSinglePoint(pos, Xvect, Yvect, Zvect, coneRM)
        drawArrow3 = @(x,y,z,varargin) quiver3(x(1),y(1),z(1),x(2)-x(1),y(2)-y(1),z(2)-z(1),0, varargin{:});
        vector_x = drawArrow3([pos(1),pos(1)+Xvect(1)],...
            [pos(2),pos(2)+Xvect(2)],...
            [pos(3),pos(3)+Xvect(3)],...
            'linewidth',2.5,'color','b');
        vector_y = drawArrow3([pos(1),pos(1)+Yvect(1)],...
            [pos(2),pos(2)+Yvect(2)],...
            [pos(3),pos(3)+Yvect(3)],...
            'linewidth',2.5,'color','r');
        vector_z = drawArrow3([pos(1),pos(1)+Zvect(1)],...
            [pos(2),pos(2)+Zvect(2)],...
            [pos(3),pos(3)+Zvect(3)],...
            'linewidth',2.5,'color','g');
        
        if plot_cone
            cone_length = 2; % in meters
            cone_angle = 5*pi()/180; % in radians
            max_radius = tan(cone_angle);
            [ax, ang] = rotMat2axisAngle(coneRM);
            
            [Xcone,Ycone,Zcone] = cylinder([0,max_radius],100); % cone of length 1
            trans=makehgtform('translate',[pos(1),pos(2),pos(3)],...
                'axisrotate', ax, ang,...
                'scale', cone_length);
            surf(Xcone,Ycone,Zcone, 'Parent', hgtransform('Matrix',trans),...
                'LineStyle', 'none', 'FaceAlpha', 0.25);
        end
    end
end