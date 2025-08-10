function R = axisAngle2rotMat(axis, angle)
% Calculate the rotation matrix associated to the axis/angle.
% Assumes right-handed reference frame such as:
%   positive angles turn anticlockwise around the directed axis
%   negative angles turn clockwise around the directed axis
% Input:
%   axis:       - Axis of rotation [x,y,z]
%   angle:      - Angle of rotation (given in radians)
%
% Output:
%   R           - matrix of rotation to apply. ex: rotated_vect = R * vect.
axis_x = axis(1)/norm(axis,2);
axis_y = axis(2)/norm(axis,2);
axis_z = axis(3)/norm(axis,2);
c = cos(angle);
s = sin(angle);
t = 1-cos(angle);

%Calculate the associated rotation matrix:
R = [c + t*axis_x^2, t*axis_x*axis_y - axis_z*s, t*axis_x*axis_z + axis_y*s;...
    t*axis_x*axis_y + axis_z*s, c + t*axis_y^2, t*axis_y*axis_z - axis_x*s;...
    t*axis_x*axis_z - axis_y*s, t*axis_y*axis_z + axis_x*s, c + t*axis_z^2];

% Round to zero the negligable elements
tolerance = 0.0001;
for i = 1:9
    if abs(R(i)) < tolerance
        R(i) = 0;
    end
end
end