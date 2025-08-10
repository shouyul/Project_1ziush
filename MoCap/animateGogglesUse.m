function animateGogglesUse(RB, latencies, opts)

if opts.record
    lapseBetweenFrames = 10; % Makes the video more fluid but the online plot lags
else
    lapseBetweenFrames = 25; % No recording so privilege visual online plot
end

hold on
%% Plot table
tabLen=1.2;
tabWid=0.8;
tabDep=0.02;
[Xtab, Ytab] = meshgrid(opts.tableLeftAngle(1):0.1:opts.tableLeftAngle(1)+tabLen,...
    opts.tableLeftAngle(2):0.1:opts.tableLeftAngle(2)+tabWid);
Ztab = ones(size(Xtab))*opts.tableLeftAngle(3);
surf(Xtab,Ytab,Ztab, 'LineStyle', 'none', 'FaceAlpha', 1, 'FaceColor', [0.94,0.87,0.80]);

%% Plot tumbler
if opts.plot_tumbler
    tumbHei=0.06;
    tumbDiam=0.06;
    [Xtumb,Ytumb,Ztumb] = cylinder(tumbDiam/2,100); % cylinder of length 1m
    trans1=makehgtform('translate',[opts.tableLeftAngle(1)+opts.tumblerWRTtable(1),...
        opts.tableLeftAngle(2)+opts.tumblerWRTtable(2),opts.tableLeftAngle(3)+opts.tumblerWRTtable(3)],...
        'scale', [1,1,tumbHei]);
    surf(Xtumb,Ytumb,Ztumb, 'Parent', hgtransform('Matrix',trans1),...
        'LineStyle', 'none', 'FaceAlpha', 1, 'FaceColor', [0.5,0.5,0.5]);
    % cap for the tumbler
    patch('XData',opts.tableLeftAngle(1)+opts.tumblerWRTtable(1)+Xtumb(2,:),...
        'YData',opts.tableLeftAngle(2)+opts.tumblerWRTtable(2)+Ytumb(2,:),...
        'ZData',opts.tableLeftAngle(3)+opts.tumblerWRTtable(3)+Ztumb(2,:)*tumbHei,...
        'LineStyle', 'none', 'FaceAlpha', 1, 'FaceColor', [0.5,0.5,0.5]);
end

%% Plot RB and cone
% Compute Rotation matrices for cone plotting
cone_RM = cat(2,reshape(RB.Zvector,size(RB.pos,1),1,size(RB.pos,2)),...
    reshape(RB.Xvector,size(RB.pos,1),1,size(RB.pos,2)),...
    reshape(RB.Yvector,size(RB.pos,1),1,size(RB.pos,2)));

% Scale vector length (default = 1m)
RB.Xvector = RB.Xvector.*0.1;
RB.Yvector = RB.Yvector.*0.25;
RB.Zvector = RB.Zvector.*0.1;

if opts.record
    % grid on produce some flickering on the video recording
    grid off
    xticks([]);
    yticks([]);
    zticks([]);
else
    grid on
    xlabel('-X');
    ylabel('Z');
    zlabel('Y');
end
axis equal

curve_RB = animatedline('Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
% Plotting loop
frameRate = (opts.animationSpeed*opts.srate)/lapseBetweenFrames;
if opts.plot_tumbler && isfield(opts,'sampStepVisData')
    latsVis = latencies(1:opts.sampStepVisData:end);
end
lats2plot = latencies(1:lapseBetweenFrames:end);
for t = 1:length(lats2plot)
    addpoints(curve_RB, RB.pos(1,lats2plot(t)), RB.pos(2,lats2plot(t)), RB.pos(3,lats2plot(t)));
    [vect_x, vect_y, vect_z, cone] = plotSinglePoint(...
        RB.pos(:,lats2plot(t)), RB.Xvector(:,lats2plot(t)), RB.Yvector(:,lats2plot(t)), RB.Zvector(:,lats2plot(t)), squeeze(cone_RM(:,:,lats2plot(t))));
    drawnow
    
    if t == 1
        textBox = uicontrol('style','text');
        set(textBox,'Units','normalized');
        set(textBox,'Position',[0.3,0.075,0.4,0.05]);
        set(textBox,'FontSize',14);
        if opts.plot_tumbler && isfield(opts,'percTumbVis')
            set(textBox,'String',...
                sprintf('t=%4.1fs: %5.1f%% visible',...
                (lats2plot(t)-lats2plot(1))/opts.srate,...
                opts.percTumbVis(find(latsVis<=lats2plot(t),1,'last'))));
        else
            set(textBox,'String',...
                sprintf('t=%4.1fs',...
                (lats2plot(t)-lats2plot(1))/opts.srate));
        end
        
        xlim([opts.tableLeftAngle(1),opts.tableLeftAngle(1)+tabLen]);
        ylim([-inf,opts.tableLeftAngle(2)+tabWid]);
        zlim([opts.tableLeftAngle(3),inf]);
        view(-40,30);
        axis vis3d
        set(gcf,'Position', [50, 50, 500, 550]);
    else
        if opts.plot_tumbler && isfield(opts,'percTumbVis')
            set(textBox,'String',...
                sprintf('t=%4.1fs: %5.1f%% visible',...
                (lats2plot(t)-lats2plot(1))/opts.srate,...
                opts.percTumbVis(find(latsVis<=lats2plot(t),1,'last'))));
        else
            set(textBox,'String',...
                sprintf('t=%4.1fs',...
                (lats2plot(t)-lats2plot(1))/opts.srate));
        end
        elapsedOneFrame = toc;
    end
    
    if t ~= length(lats2plot)
        if t == 1
            pause(1/frameRate);
        else
            if elapsedOneFrame <= 1/frameRate
                pause(1/frameRate - elapsedOneFrame);
            else
                %fprintf('Refreshing rate too high (%.3fs): last frame took %.3fs to display.\n',...
                %    1/frameRate, elapsedOneFrame)
            end
        end
        tic;
    end
    
    if opts.record
        Frames(t) = getframe(gcf);
    end
    
    if t ~= length(lats2plot)
        delete(vect_x);
        delete(vect_y);
        delete(vect_z);
        delete(cone);
    end
end

if opts.record
    video = VideoWriter(opts.titleFile, 'MPEG-4');
    video.Quality = 100;
    video.FrameRate = frameRate;
    %video.FrameRate = (length(Frames)-1)*opts.srate/(latencies(end) - latencies(1));
    open(video)
    writeVideo(video, Frames);
    close(video)
    
    close gcf
end

    function [vector_x, vector_y, vector_z, cone] = plotSinglePoint(pos, Xvect, Yvect, Zvect, coneRM)
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
        
        cone_length = 2; % in meters
        cone_angle = 5*pi()/180; % in radians
        max_radius = tan(cone_angle);
        [ax, ang] = rotMat2axisAngle(coneRM);
        
        [Xcone,Ycone,Zcone] = cylinder([0,max_radius],100); % cone of length 1
        trans=makehgtform('translate',[pos(1),pos(2),pos(3)],...
            'axisrotate', ax, ang,...
            'scale', cone_length);
        cone = surf(Xcone,Ycone,Zcone, 'Parent', hgtransform('Matrix',trans),...
            'LineStyle', 'none', 'FaceAlpha', 0.25);
    end
end