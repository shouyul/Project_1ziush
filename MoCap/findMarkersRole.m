function [EEG_mocap2, badSamples] = findMarkersRole(EEG_mocap, rb, Nmarkers, tol_params)
EEG_mocap2 = EEG_mocap;
badSamples = false(1,EEG_mocap.pnts);

% Find blocks to inspect for each rb
switch rb
    case 'RB1'
        bl2inspect = unique(EEG_mocap.etc.TrialData.Block(:));
        accepted_roles = {'FrontRight','TopRight','FrontLeft','BackLeft'};
    case 'RB2'
        bl2inspect = unique(EEG_mocap.etc.TrialData.Block(strcmp(EEG_mocap.etc.TrialData.Condition, 'GogglesON')));
        accepted_roles = {'Right','Left','Bottom','Camera'};
    otherwise
        error('Unknown rigid body')
end

% Some useful variables
rb_chans = contains({EEG_mocap.chanlocs.labels},rb);
%bnd_evts = find(strcmp({EEG_mocap.event.type}, 'boundary'));

%% Compute latencies to inspect
all_latencies = 1:EEG_mocap2.pnts;
latencies2inspect = all_latencies(~all(EEG_mocap2.data(rb_chans,:)==0, 1));
n_pnts2inspect = length(latencies2inspect);

%%%%%%%% First assignment %%%%%%%%%
%% Plot the markers relative position
plotMarkers(EEG_mocap, rb_chans, Nmarkers, latencies2inspect(1));

%% Ask for new names to the user
fprintf('Accepted roles: ');
for a = 1:numel(accepted_roles)
    if a < numel(accepted_roles)
        fprintf('%s, ',accepted_roles{a})
    else
        fprintf('%s\n',accepted_roles{a})
    end
end

newNames = assignRolesLoop(accepted_roles);
%% Reassign names in the chanlocs structure
for m = 1:Nmarkers
    Mk_chans = contains({EEG_mocap.chanlocs.labels},sprintf('Marker%d',m-1));
    chans2change = find(rb_chans&Mk_chans);
    for c = chans2change
        EEG_mocap2.chanlocs(c).labels = replace(EEG_mocap.chanlocs(c).labels,...
            sprintf('Marker%d',m-1), newNames.(sprintf('Marker%d',m-1)));
    end
end

%% Check validity and compute Rigid body to understand markers relative position
% Now check dimensions of the headset (if RB2)
if ~strcmp(rb, 'RB2') || checkDimensionsRB2(EEG_mocap2, rb_chans, latencies2inspect(1), tol_params.dims)
    % Sample validated, continue
    [~, RB_ref] = rigidBodyProperties(EEG_mocap2, rb, latencies2inspect(1));
    forcedInspection = false;
else
    % Sample rejected
    fprintf('Latency %d: dimensions violated.\n', latencies2inspect(1));
    [~, RB_ref_bad{1}] = rigidBodyProperties(EEG_mocap2, rb, latencies2inspect(1));
    RB_ref_bad_last = RB_ref_bad{1};
    badSamples(latencies2inspect(1)) = true;
    forcedInspection = true;
end

progressBar = waitbar(1/n_pnts2inspect, 'Checking markers role point by point',...
    'Name', sprintf('%s - findMarkersRole',rb));
for p = 2:n_pnts2inspect
    lat = latencies2inspect(p);
    
    if lat ~= latencies2inspect(p-1)+1
        % Break in latencies, force new inspection
        forcedInspection = true;
    end
    
    if ~forcedInspection
        % Consecutive latencies
        [~, RB] = rigidBodyProperties(EEG_mocap2, rb, lat);
        
        % Try to assign markers automatically by correspondance to the reference RB
        [corrDirMat, corrDistMat, assignment] = autoAssignment(RB, RB_ref, newNames, tol_params.ori_good, tol_params.dist_good);
        
        if all(diag(assignment))
            % The marker roles seem respected
            % Now check dimensions of the headset (if RB2)
            if ~strcmp(rb, 'RB2') || checkDimensionsRB2(EEG_mocap2, rb_chans, lat, tol_params.dims)
                % Sample validated
                RB_ref = RB; % Change RB_ref to keep as close as possible to the current spatial organization
            else
                % Sample rejected
                badSamples(lat) = true;
                
                % Automatic check to see if it looks like an already registered bad sample configuration
                % (to avoid oversizing RB_ref_bad to save computation time)
                if exist('RB_ref_bad','var')
                    found_bad_ref = false;
                    % First try with RB_ref_bad_last (may save computation time)
                    [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad_last,...
                        newNames, tol_params.ori_bad, tol_params.dist_bad);
                    if all(diag(assignment_bad))
                        fprintf('Latency %d: bad sample detected automatically (last ref).\n',lat)
                        found_bad_ref = true;
                    end
                    
                    if ~found_bad_ref
                        % Then loop over all bad references
                        for b = numel(RB_ref_bad):-1:1
                            [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad{b},...
                                newNames, tol_params.ori_bad, tol_params.dist_bad);
                            
                            if all(diag(assignment_bad))
                                fprintf('Latency %d: bad sample detected automatically (ref #%d).\n',lat,b)
                                RB_ref_bad_last = RB_ref_bad{b};
                                found_bad_ref = true;
                                break
                            end
                        end
                    end
                    
                    if ~found_bad_ref
                        fprintf('Latency %d: dimensions violated.\n', lat);
                        % Assign new Bad reference
                        RB_ref_bad{numel(RB_ref_bad)+1} = RB;
                        RB_ref_bad_last = RB;
                    end
                else
                    fprintf('Latency %d: dimensions violated.\n', lat);
                    RB_ref_bad{1} = RB;
                    RB_ref_bad_last = RB;
                end
            end
        else
            if all(sum(assignment,2)==1)
                fprintf('Latency %d: roles swap detected.\n', lat);
                % No confusion, reassign
                for m1 = 1:Nmarkers
                    %if maxCorresp(m1) == m1
                    if assignment(m1,m1)
                        % This marker kept its role, no need to exchange data
                    else
                        Mk_chans = contains({EEG_mocap.chanlocs.labels}, sprintf('Marker%d',m1-1));
                        %Mk_chans2 = contains({EEG_mocap2.chanlocs.labels}, newNames.(sprintf('Marker%d',maxCorresp(m1)-1)));
                        Mk_chans2 = contains({EEG_mocap2.chanlocs.labels}, newNames.(sprintf('Marker%d',find(assignment(m1,:))-1)));
                        % Replace data in EEG_mocap2
                        EEG_mocap2.data(rb_chans&Mk_chans2,lat:end) = EEG_mocap.data(rb_chans&Mk_chans,lat:end);
                    end
                end
                % Replace data in EEG_mocap to copy the right data later
                EEG_mocap.data(rb_chans,lat:end) = EEG_mocap2.data(rb_chans,lat:end);
                
                % Now check dimensions of the headset (if RB2)
                if ~strcmp(rb, 'RB2') || checkDimensionsRB2(EEG_mocap2, rb_chans, lat, tol_params.dims)
                    % Sample validated
                    % Recompute ref
                    [~, RB_ref] = rigidBodyProperties(EEG_mocap2, rb, lat);
                else
                    % Sample rejected
                    badSamples(lat) = true;
                    
                    % Automatic check to see if it looks like an already registered bad sample configuration
                    % (to avoid oversizing RB_ref_bad to save computation time)
                    [~, RB] = rigidBodyProperties(EEG_mocap2, rb, lat);
                    if exist('RB_ref_bad','var')
                        found_bad_ref = false;
                        % First try with RB_ref_bad_last (may save computation time)
                        [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad_last,...
                            newNames, tol_params.ori_bad, tol_params.dist_bad);
                        if all(diag(assignment_bad))
                            fprintf('Latency %d: bad sample detected automatically (last ref).\n',lat)
                            found_bad_ref = true;
                        end
                        
                        if ~found_bad_ref
                            % Then loop over all bad references
                            for b = numel(RB_ref_bad):-1:1
                                [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad{b},...
                                    newNames, tol_params.ori_bad, tol_params.dist_bad);
                                
                                if all(diag(assignment_bad))
                                    fprintf('Latency %d: bad sample detected automatically (ref #%d).\n',lat,b)
                                    RB_ref_bad_last = RB_ref_bad{b};
                                    found_bad_ref = true;
                                    break
                                end
                            end
                        end
                        
                        if ~found_bad_ref
                            fprintf('Latency %d: dimensions violated.\n', lat);
                            % Assign new Bad reference
                            RB_ref_bad{numel(RB_ref_bad)+1} = RB;
                            RB_ref_bad_last = RB;
                        end
                    else
                        fprintf('Latency %d: dimensions violated.\n', lat);
                        RB_ref_bad{1} = RB;
                        RB_ref_bad_last = RB;
                    end
                end
            else
                fprintf('Latency %d: Could not assign all roles automatically.\n', lat);
                %% Automatic check to see if it looks like a bad sample configuration
                if exist('RB_ref_bad','var')
                    found_bad_ref = false;
                    % First try with RB_ref_bad_last (may save computation time)
                    [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad_last,...
                        newNames, tol_params.ori_bad, tol_params.dist_bad);
                    if all(diag(assignment_bad))
                        fprintf('Latency %d: bad sample detected automatically (last ref).\n',lat)
                        found_bad_ref = true;
                        badSamples(lat) = true;
                    end
                    
                    if ~found_bad_ref
                        % Then loop over all bad references
                        for b = numel(RB_ref_bad):-1:1
                            [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad{b},...
                                newNames, tol_params.ori_bad, tol_params.dist_bad);
                            
                            if all(diag(assignment_bad))
                                fprintf('Latency %d: bad sample detected automatically (ref #%d).\n',lat,b)
                                RB_ref_bad_last = RB_ref_bad{b};
                                badSamples(lat) = true;
                                break
                            end
                        end
                    end
                end
                
                if ~badSamples(lat)
                    % Not labelled bad until now
                    %% Ask user
                    plotMarkers(EEG_mocap, rb_chans, Nmarkers, lat);
                    fprintf('CorrDirMat                 CorrDistMat              Assignment\n')
                    for j = 1:Nmarkers
                        fprintf('%5.2f %5.2f %5.2f %5.2f   %5.2f %5.2f %5.2f %5.2f   %d %d %d %d\n',...
                            corrDirMat(j,1), corrDirMat(j,2), corrDirMat(j,3), corrDirMat(j,4),...
                            corrDistMat(j,1), corrDistMat(j,2), corrDistMat(j,3), corrDistMat(j,4),...
                            assignment(j,1), assignment(j,2), assignment(j,3), assignment(j,4))
                    end
                    valid = input('Is it a valid position? [0/1] ');
                    if valid
                        change = input('Change roles? [0/1] ');
                        if change
                            newNames2 = assignRolesLoop(accepted_roles);
                            %% Reassign
                            for m = 1:Nmarkers
                                if ~strcmp(newNames.(sprintf('Marker%d',m-1)), newNames2.(sprintf('Marker%d',m-1)))
                                    % Find channels correspondance
                                    Mk_chans = contains({EEG_mocap.chanlocs.labels}, sprintf('Marker%d',m-1));
                                    Mk_chans2 = contains({EEG_mocap2.chanlocs.labels}, newNames2.(sprintf('Marker%d',m-1)));
                                    
                                    % Replace data in EEG_mocap2
                                    EEG_mocap2.data(rb_chans&Mk_chans2,lat:end) = EEG_mocap.data(rb_chans&Mk_chans,lat:end);
                                end
                            end
                            % Replace data in EEG_mocap to copy the right data later
                            EEG_mocap.data(rb_chans,lat:end) = EEG_mocap2.data(rb_chans,lat:end);
                            
                            % Now check dimensions of the headset (if RB2)
                            if ~strcmp(rb, 'RB2') || checkDimensionsRB2(EEG_mocap2, rb_chans, lat, tol_params.dims)
                                % Recompute ref
                                [~, RB_ref] = rigidBodyProperties(EEG_mocap2, rb, lat);
                                close gcf
                            else
                                % Sample rejected
                                badSamples(lat) = true;
                                
                                % Automatic check to see if it looks like an already registered bad sample configuration
                                % (to avoid oversizing RB_ref_bad to save computation time)
                                [~, RB] = rigidBodyProperties(EEG_mocap2, rb, lat);
                                if exist('RB_ref_bad','var')
                                    found_bad_ref = false;
                                    % First try with RB_ref_bad_last (may save computation time)
                                    [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad_last,...
                                        newNames, tol_params.ori_bad, tol_params.dist_bad);
                                    if all(diag(assignment_bad))
                                        fprintf('Latency %d: bad sample detected automatically (last ref).\n',lat)
                                        found_bad_ref = true;
                                    end
                                    
                                    if ~found_bad_ref
                                        % Then loop over all bad references
                                        for b = numel(RB_ref_bad):-1:1
                                            [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad{b},...
                                                newNames, tol_params.ori_bad, tol_params.dist_bad);
                                            
                                            if all(diag(assignment_bad))
                                                fprintf('Latency %d: bad sample detected automatically (ref #%d).\n',lat,b)
                                                RB_ref_bad_last = RB_ref_bad{b};
                                                found_bad_ref = true;
                                                close gcf
                                                break
                                            end
                                        end
                                    end
                                    
                                    if ~found_bad_ref
                                        fprintf('Latency %d: dimensions violated.\n', lat);
                                        % Assign new Bad reference
                                        RB_ref_bad{numel(RB_ref_bad)+1} = RB;
                                        RB_ref_bad_last = RB;
                                        title(gca, sprintf('Latency %d - Bad reference #%d',lat, numel(RB_ref_bad)));
                                    end
                                else
                                    fprintf('Latency %d: dimensions violated.\n', lat);
                                    RB_ref_bad{1} = RB;
                                    RB_ref_bad_last = RB;
                                    title(gca, sprintf('Latency %d - Bad reference #%d',lat, numel(RB_ref_bad)));
                                end
                            end
                        else
                            % Now check dimensions of the headset (if RB2)
                            if ~strcmp(rb, 'RB2') || checkDimensionsRB2(EEG_mocap2, rb_chans, lat, tol_params.dims)
                                RB_ref = RB; % Change RB_ref to keep as close as possible to the current spatial organization
                                close gcf
                            else
                                % Sample rejected
                                badSamples(lat) = true;
                                
                                % No need to check bad refs (already done earlier and no change of markers role)
                                fprintf('Latency %d: dimensions violated.\n', lat);
                                % Assign new Bad reference
                                if ~exist('RB_ref_bad','var')
                                    RB_ref_bad{1} = RB;
                                    RB_ref_bad_last = RB;
                                else
                                    RB_ref_bad{numel(RB_ref_bad)+1} = RB;
                                    RB_ref_bad_last = RB;
                                end
                                title(gca, sprintf('Latency %d - Bad reference #%d',lat, numel(RB_ref_bad)));
                            end
                        end
                    else
                        fprintf('Latency %d: labelled bad by user.\n', lat);
                        badSamples(lat) = true;
                        % Assign new Bad reference
                        if ~exist('RB_ref_bad','var')
                            RB_ref_bad{1} = RB;
                            RB_ref_bad_last = RB;
                        else
                            RB_ref_bad{numel(RB_ref_bad)+1} = RB;
                            RB_ref_bad_last = RB;
                        end
                        title(gca, sprintf('Latency %d - Bad reference #%d',lat, numel(RB_ref_bad)));
                    end
                end
            end
        end
    else
        %% Forced Inspection
        plotMarkers(EEG_mocap, rb_chans, Nmarkers, lat);
        valid = input('Is it a valid position? [0/1] ');
        if valid
            newNames2 = assignRolesLoop(accepted_roles);
            % Reassign
            for m = 1:Nmarkers
                if ~strcmp(newNames.(sprintf('Marker%d',m-1)), newNames2.(sprintf('Marker%d',m-1)))
                    % Find channels correspondance
                    Mk_chans = contains({EEG_mocap.chanlocs.labels}, sprintf('Marker%d',m-1));
                    Mk_chans2 = contains({EEG_mocap2.chanlocs.labels}, newNames2.(sprintf('Marker%d',m-1)));
                    
                    % Replace data in EEG_mocap2
                    EEG_mocap2.data(rb_chans&Mk_chans2,lat:end) = EEG_mocap.data(rb_chans&Mk_chans,lat:end);
                end
            end
            % Replace data in EEG_mocap to copy the right data later
            EEG_mocap.data(rb_chans,lat:end) = EEG_mocap2.data(rb_chans,lat:end);
            
            % Now check dimensions of the headset (if RB2)
            if ~strcmp(rb, 'RB2') || checkDimensionsRB2(EEG_mocap2, rb_chans, lat, tol_params.dims)
                % Sample validated, continue
                [~, RB_ref] = rigidBodyProperties(EEG_mocap2, rb, lat);
                forcedInspection = false;
                close gcf
            else
                badSamples(lat) = true;               
                
                % Automatic check to see if it looks like an already registered bad sample configuration
                % (to avoid oversizing RB_ref_bad to save computation time)
                [~, RB] = rigidBodyProperties(EEG_mocap2, rb, lat);
                if exist('RB_ref_bad','var')
                    found_bad_ref = false;
                    % First try with RB_ref_bad_last (may save computation time)
                    [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad_last,...
                        newNames, tol_params.ori_bad, tol_params.dist_bad);
                    if all(diag(assignment_bad))
                        fprintf('Latency %d: bad sample detected automatically (last ref).\n',lat)
                        found_bad_ref = true;
                    end
                    
                    if ~found_bad_ref
                        % Then loop over all bad references
                        for b = numel(RB_ref_bad):-1:1
                            [~, ~, assignment_bad] = autoAssignment(RB, RB_ref_bad{b},...
                                newNames, tol_params.ori_bad, tol_params.dist_bad);
                            
                            if all(diag(assignment_bad))
                                fprintf('Latency %d: bad sample detected automatically (ref #%d).\n',lat,b)
                                RB_ref_bad_last = RB_ref_bad{b};
                                found_bad_ref = true;
                                close gcf
                                break
                            end
                        end
                    end
                    
                    if ~found_bad_ref
                        fprintf('Latency %d: dimensions violated.\n', lat);
                        % Assign new Bad reference
                        RB_ref_bad{numel(RB_ref_bad)+1} = RB;
                        RB_ref_bad_last = RB;
                        title(gca, sprintf('Latency %d - Bad reference #%d',lat, numel(RB_ref_bad)));
                    end
                else
                    fprintf('Latency %d: dimensions violated.\n', lat);
                    RB_ref_bad{1} = RB;
                    RB_ref_bad_last = RB;
                    title(gca, sprintf('Latency %d - Bad reference #%d',lat, numel(RB_ref_bad)));
                end
            end
        else
            fprintf('Latency %d: labelled bad by user.\n', lat);
            badSamples(lat) = true;
            % Assign new Bad reference
            if ~exist('RB_ref_bad','var')
                [~, RB_ref_bad{1}] = rigidBodyProperties(EEG_mocap2, rb, lat);
                RB_ref_bad_last = RB_ref_bad{1};
            else
                [~, RB_ref_bad{numel(RB_ref_bad)+1}] = rigidBodyProperties(EEG_mocap2, rb, lat);
                RB_ref_bad_last = RB_ref_bad{numel(RB_ref_bad)};
            end
            title(gca, sprintf('Latency %d - Bad reference #%d',lat, numel(RB_ref_bad)));
        end
    end
    waitbar(p/n_pnts2inspect, progressBar,...
        sprintf('Checking markers role point by point: %.1f%%',100*p/n_pnts2inspect));
end
delete(progressBar);



% firstBlock = true;
% for bl = bl2inspect'
%     blStart = intersect(find(strcmp({EEG_mocap.event.type}, 'BlockStart')),...
%         findInStructWithEmpties(EEG_mocap.event, 'block', bl));
%     latencies = ceil(EEG_mocap.event(blStart).latency-EEG_mocap.srate/2):...
%         floor(EEG_mocap.event(blStart).latency+EEG_mocap.srate/2);
%
%     %% Plot the markers relative position
%     figure;
%     hold on
%     leg = [];
%     leg_labels = {};
%     for m = 1:Nmarkers
%         % Plot the average position of markers at the beginning of the recording
%         Mk_chans = contains({EEG_mocap.chanlocs.labels},sprintf('Marker%d',m-1));
%         Mk_pos = mean(EEG_mocap.data(rb_chans&Mk_chans,latencies),2);
%         s = scatter3(Mk_pos(1), Mk_pos(2), Mk_pos(3));
%         leg = [leg,s];
%         leg_labels = [leg_labels, sprintf('Marker%d',m-1)];
%
%         % Plot line between current marker and previous ones
%         for m_prev = 1:m-1
%             Mk_prev_chans = contains({EEG_mocap.chanlocs.labels},sprintf('Marker%d',m_prev-1));
%             Mk_prev_pos = mean(EEG_mocap.data(rb_chans&Mk_prev_chans,latencies),2);
%             plot3([Mk_prev_pos(1);Mk_pos(1)], [Mk_prev_pos(2);Mk_pos(2)],[Mk_prev_pos(3);Mk_pos(3)], 'k')
%         end
%     end
%
%     legend(leg, leg_labels);
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     axis equal
%     title(sprintf('Block %d',bl))
%
%     %% Ask for new names to the user
%     if firstBlock
%         fprintf('Accepted roles: ');
%         for a = 1:numel(accepted_roles)
%             if a < numel(accepted_roles)
%                 fprintf('%s, ',accepted_roles{a})
%             else
%                 fprintf('%s\n',accepted_roles{a})
%             end
%         end
%
%         newNames = assignRolesLoop(accepted_roles);
%
%         %% Reassign names in the chanlocs structure
%         for m = 1:Nmarkers
%             Mk_chans = contains({EEG_mocap.chanlocs.labels},sprintf('Marker%d',m-1));
%             chans2change = find(rb_chans&Mk_chans);
%             for c = chans2change
%                 EEG_mocap2.chanlocs(c).labels = replace(EEG_mocap.chanlocs(c).labels,...
%                     sprintf('Marker%d',m-1), newNames.(sprintf('Marker%d',m-1)));
%             end
%         end
%         firstBlock = false;
%     else
%         fprintf('Assigned roles:\n');
%         for m = 1:Nmarkers
%             fprintf('%s: %s\n',sprintf('Marker%d',m-1),newNames.(sprintf('Marker%d',m-1)))
%         end
%         change = input('Did they change for this block? [0/1] ');
%         if change
%             newNames2 = assignRolesLoop(accepted_roles);
%
%             for m = 1:Nmarkers
%                 if ~strcmp(newNames.(sprintf('Marker%d',m-1)), newNames2.(sprintf('Marker%d',m-1)))
%                     % Find channels correspondance
%                     old_Mk_chans = contains({EEG_mocap.chanlocs.labels}, sprintf('Marker%d',m-1));
%                     new_Mk_chans = contains({EEG_mocap2.chanlocs.labels},...
%                         sprintf('_%s_', newNames2.(sprintf('Marker%d',m-1)))); % may not be unique because of contains
%
%                     % Find latencies to invert
%                     start_ev = find([EEG_mocap.event(bnd_evts).latency] < EEG_mocap.event(blStart).latency, 1, 'last');
%                     if ~isempty(start_ev)
%                         start_bl = ceil(EEG_mocap.event(bnd_evts(start_ev)).latency);
%                     else
%                         start_bl = 1;
%                     end
%                     stop_ev = find([EEG_mocap.event(bnd_evts).latency] > EEG_mocap.event(blStart).latency, 1);
%                     if ~isempty(stop_ev)
%                         stop_bl = floor(EEG_mocap.event(bnd_evts(stop_ev)).latency);
%                     else
%                         stop_bl = EEG_mocap.pnts;
%                     end
%                     lats_bl = start_bl:stop_bl;
%
%                     % Replace data in EEG_mocap2
%                     EEG_mocap2.data(rb_chans&new_Mk_chans&X_chans,lats_bl) =...
%                         EEG_mocap.data(rb_chans&old_Mk_chans&X_chans,lats_bl);
%                     EEG_mocap2.data(rb_chans&new_Mk_chans&Y_chans,lats_bl) =...
%                         EEG_mocap.data(rb_chans&old_Mk_chans&Y_chans,lats_bl);
%                     EEG_mocap2.data(rb_chans&new_Mk_chans&Z_chans,lats_bl) =...
%                         EEG_mocap.data(rb_chans&old_Mk_chans&Z_chans,lats_bl);
%                 end
%             end
%         else
%             % nothing to do
%         end
%     end
%end

    function plotMarkers(EEG_mocap, rb_chans, Nmarkers, lats)
        % lats can be single lat or mutliple lats (plotting the average in that case)
        X_chs = contains({EEG_mocap.chanlocs.labels},'X');
        Y_chs = contains({EEG_mocap.chanlocs.labels},'Y');
        Z_chs = contains({EEG_mocap.chanlocs.labels},'Z');
        
        figure;
        hold on
        leg = [];
        leg_labels = {};
        for i = 1:Nmarkers
            % Plot the average position of markers at the beginning of the recording
            Mk_chs = contains({EEG_mocap.chanlocs.labels},sprintf('Marker%d',i-1));
            Mk_pos = mean(EEG_mocap.data(rb_chans&Mk_chs,lats),2);
            s = scatter3(-Mk_pos(1), Mk_pos(3), Mk_pos(2));
            leg = [leg,s];
            leg_labels = [leg_labels, sprintf('Marker%d',i-1)];
            
            % Plot line between current marker and previous ones
            for m_prev = 1:i-1
                Mk_prev_chans = contains({EEG_mocap.chanlocs.labels},sprintf('Marker%d',m_prev-1));
                Mk_prev_pos = mean(EEG_mocap.data(rb_chans&Mk_prev_chans,lats),2);
                plot3(-[Mk_prev_pos(1);Mk_pos(1)], [Mk_prev_pos(3);Mk_pos(3)],[Mk_prev_pos(2);Mk_pos(2)], 'k')
            end
        end
        
        %% Additionally plot centroid (better 3D vizualization)
        centroid = [mean(mean(EEG_mocap.data(rb_chans&X_chs,lats),1),2);...
            mean(mean(EEG_mocap.data(rb_chans&Y_chs,lats),1),2);...
            mean(mean(EEG_mocap.data(rb_chans&Z_chs,lats),1),2)];
        
        s = scatter3(-centroid(1), centroid(3), centroid(2), 'k');
        leg = [leg,s];
        leg_labels = [leg_labels, 'Centroid'];
        
        % Plot line between current marker and previous ones
        for i = 1:Nmarkers
            Mk_chs = contains({EEG_mocap.chanlocs.labels},sprintf('Marker%d',i-1));
            Mk_pos = mean(EEG_mocap.data(rb_chans&Mk_chs,lats),2);
            plot3(-[centroid(1);Mk_pos(1)], [centroid(3);Mk_pos(3)],[centroid(2);Mk_pos(2)], '--k')
        end
        
        legend(leg, leg_labels);
        xlabel('-X');
        ylabel('Z');
        zlabel('Y');
        axis equal
        xlim([0,1.2]);
        ylim([-inf,0.2]);
        zlim([0,inf]);
        view(-125,15);
        grid on
        title(sprintf('Latency %d', mean(lats)));
    end

    function newNames = assignRolesLoop(accepted_roles)
        newNames = struct();
        remaining_roles = accepted_roles;
        for i = 1:numel(accepted_roles)
            fprintf('Remaining roles: ');
            for r = 1:numel(remaining_roles)
                if r < numel(remaining_roles)
                    fprintf('%s, ',remaining_roles{r})
                else
                    fprintf('%s\n',remaining_roles{r})
                end
            end
            
            while true
                role = input(sprintf('Role of marker %d? ', i-1));
                if ischar(role)
                    if any(strcmp(remaining_roles, role))
                        newNames.(sprintf('Marker%d',i-1)) = role;
                        remaining_roles = remaining_roles(~strcmp(remaining_roles,role));
                        break
                    else
                        fprintf('%s not understood\n', role)
                        continue
                    end
                else
                    disp('Role should be a char')
                    continue
                end
            end
        end
    end

    function [centroid, RB] = rigidBodyProperties(EEG_mocap, rb, lats)
        % check inputs
        if nargin < 3 || (ischar(lats) && stcrmp(lats, 'all'))
            lats = 1:EEG_mocap.pnts;
        end
        
        rb_chs = contains({EEG_mocap.chanlocs.labels}, rb);
        X_chs = contains({EEG_mocap.chanlocs.labels}, 'X');
        Y_chs = contains({EEG_mocap.chanlocs.labels}, 'Y');
        Z_chs = contains({EEG_mocap.chanlocs.labels}, 'Z');
        
        centroid = [mean(EEG_mocap.data(rb_chs&X_chs,lats),1);...
            mean(EEG_mocap.data(rb_chs&Y_chs,lats),1);...
            mean(EEG_mocap.data(rb_chs&Z_chs,lats),1)];
        
        switch rb
            case 'RB1'
                FR_chs = contains({EEG_mocap.chanlocs.labels}, 'FrontRight');
                RB.FrontRight = EEG_mocap.data(rb_chs&FR_chs,lats) - centroid;
                TR_chs = contains({EEG_mocap.chanlocs.labels}, 'TopRight');
                RB.TopRight = EEG_mocap.data(rb_chs&TR_chs,lats) - centroid;
                FL_chs = contains({EEG_mocap.chanlocs.labels}, 'FrontLeft');
                RB.FrontLeft = EEG_mocap.data(rb_chs&FL_chs,lats) - centroid;
                BL_chs = contains({EEG_mocap.chanlocs.labels}, 'BackLeft');
                RB.BackLeft = EEG_mocap.data(rb_chs&BL_chs,lats) - centroid;
            case 'RB2'
                R_chs = contains({EEG_mocap.chanlocs.labels}, 'Right');
                RB.Right = EEG_mocap.data(rb_chs&R_chs,lats) - centroid;
                L_chs = contains({EEG_mocap.chanlocs.labels}, 'Left');
                RB.Left = EEG_mocap.data(rb_chs&L_chs,lats) - centroid;
                C_chs = contains({EEG_mocap.chanlocs.labels}, 'Camera');
                RB.Camera = EEG_mocap.data(rb_chs&C_chs,lats) - centroid;
                B_chs = contains({EEG_mocap.chanlocs.labels}, 'Bottom');
                RB.Bottom = EEG_mocap.data(rb_chs&B_chs,lats) - centroid;
        end
    end

    function [corrDirMat, corrDistMat, assignment] = autoAssignment(RB, RB_ref, newNames, tol_ori, tol_dist)
        Nmark = numel(fieldnames(newNames));
        corrDirMat = zeros(Nmark);
        corrDistMat = ones(Nmark);
        for mk1 = 1:Nmark
            for mk2 = 1:Nmark
                corrDirMat(mk1,mk2) = dot(normc(RB.(newNames.(sprintf('Marker%d',mk2-1)))),...
                    normc(RB_ref.(newNames.(sprintf('Marker%d',mk1-1)))));
                corrDistMat(mk1,mk2) = abs(norm(RB.(newNames.(sprintf('Marker%d',mk2-1)))) - ...
                    norm(RB_ref.(newNames.(sprintf('Marker%d',mk1-1)))))/norm(RB_ref.(newNames.(sprintf('Marker%d',mk1-1))));
            end
        end
        assignment = corrDirMat >= cos(tol_ori*pi()/180) & corrDistMat <= tol_dist;
    end

    function good = checkDimensionsRB2(EEG_mocap, rb_chs, lat, tol)
        distLR_ref = 0.17; %in meters
        distCB_ref = 0.056; %in meters
        
        LeftMk_chans = contains({EEG_mocap.chanlocs.labels}, 'Left');
        RightMk_chans = contains({EEG_mocap.chanlocs.labels}, 'Right');
        distLR = norm(EEG_mocap.data(rb_chs&LeftMk_chans,lat)-EEG_mocap.data(rb_chs&RightMk_chans,lat));
        CamMk_chans = contains({EEG_mocap.chanlocs.labels}, 'Camera');
        BotMk_chans = contains({EEG_mocap.chanlocs.labels}, 'Bottom');
        distCB = norm(EEG_mocap.data(rb_chs&CamMk_chans,lat)-EEG_mocap.data(rb_chs&BotMk_chans,lat));
        
        if distLR >= distLR_ref*(1-tol) && distLR <= distLR_ref*(1+tol) &&...
                distCB >= distCB_ref*(1-tol) && distCB <= distCB_ref*(1+tol)
            % Real distances are respected
            good = true;
        else
            good = false;
        end
    end
end