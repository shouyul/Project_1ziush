function RB = computeRigidBody(EEG_mocap, rb, lats)
% check inputs
if nargin < 3 || (ischar(lats) && stcrmp(lats, 'all'))
    lats = 1:EEG_mocap.pnts;
end

rb_chans = contains({EEG_mocap.chanlocs.labels}, rb);
X_chans = contains({EEG_mocap.chanlocs.labels}, 'X');
Y_chans = contains({EEG_mocap.chanlocs.labels}, 'Y');
Z_chans = contains({EEG_mocap.chanlocs.labels}, 'Z');

switch rb
    case 'RB1'
                FR_chans = contains({EEG_mocap.chanlocs.labels}, 'FrontRight');
        TR_chans = contains({EEG_mocap.chanlocs.labels}, 'TopRight');
        FL_chans = contains({EEG_mocap.chanlocs.labels}, 'FrontLeft');
        BL_chans = contains({EEG_mocap.chanlocs.labels}, 'BackLeft');
        
        % RB position as grand average of all markers
        RB.pos = [mean(EEG_mocap.data(rb_chans&X_chans,lats),1);...
            mean(EEG_mocap.data(rb_chans&Y_chans,lats),1);...
            mean(EEG_mocap.data(rb_chans&Z_chans,lats),1)];
        
        % Orientation Vectors
        RB.Xvector = normc(EEG_mocap.data(rb_chans&FR_chans,lats) - EEG_mocap.data(rb_chans&FL_chans,lats));
        vect1 = EEG_mocap.data(rb_chans&TR_chans,lats) - EEG_mocap.data(rb_chans&FL_chans,lats);
        RB.Yvector = normc(dot(RB.Xvector,vect1).*RB.Xvector - vect1);
        RB.Zvector = cross(RB.Xvector,RB.Yvector);
    case 'RB2'
        R_chans = contains({EEG_mocap.chanlocs.labels}, 'Right');
        L_chans = contains({EEG_mocap.chanlocs.labels}, 'Left');
        C_chans = contains({EEG_mocap.chanlocs.labels}, 'Camera');
        B_chans = contains({EEG_mocap.chanlocs.labels}, 'Bottom');
        
        % RB position centered on camera marker
        RB.pos = [EEG_mocap.data(rb_chans&X_chans&C_chans,lats);...
            EEG_mocap.data(rb_chans&Y_chans&C_chans,lats);...
            EEG_mocap.data(rb_chans&Z_chans&C_chans,lats)];       
        
        % Orientation Vectors
        RB.Zvector = normc(EEG_mocap.data(rb_chans&C_chans,lats) - EEG_mocap.data(rb_chans&B_chans,lats));
        vect1 = normc(EEG_mocap.data(rb_chans&R_chans,lats) - EEG_mocap.data(rb_chans&L_chans,lats));
        RB.Yvector = cross(RB.Zvector,vect1);
        RB.Xvector = cross(RB.Yvector,RB.Zvector);
        
        %         RB.Xvector = normc(EEG_mocap.data(rb_chans&R_chans,:) - EEG_mocap.data(rb_chans&L_chans,:));
        %         vect1 = EEG_mocap.data(rb_chans&C_chans,:) - EEG_mocap.data(rb_chans&L_chans,:);
        %         RB.Yvector = normc(vect1 - dot(RB.Xvector,vect1).*RB.Xvector);
        %         RB.Zvector = cross(RB.Xvector,RB.Yvector);
        
        % Correct RB.pos for marker dimensions (get closer to actual camera position)
        marker_rad = 0.007; % in meters
        RB.pos = RB.pos - marker_rad*RB.Yvector - 2*marker_rad*RB.Zvector;
end

% Reference frame conversion from Optitrack space to standard space
pos_copy = RB.pos;
Xvect_copy = RB.Xvector;
Yvect_copy = RB.Yvector;
Zvect_copy = RB.Zvector;

RB.pos = cat(1,-pos_copy(1,:),pos_copy(3,:),pos_copy(2,:));
RB.Xvector = cat(1,-Xvect_copy(1,:),Xvect_copy(3,:),Xvect_copy(2,:));
RB.Yvector = cat(1,-Yvect_copy(1,:),Yvect_copy(3,:),Yvect_copy(2,:));
RB.Zvector = cat(1,-Zvect_copy(1,:),Zvect_copy(3,:),Zvect_copy(2,:));
end