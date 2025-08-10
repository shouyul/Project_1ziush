function [axis, angle] = rotMat2axisAngle(R)
% Calculate the axis/angle associated to the rotation matrix.
% Assumes right-handed reference frame such as:
%   positive angles turn anticlockwise around the directed axis
%   negative angles turn clockwise around the directed axis
% Input:
%   R           - matrix of rotation to apply. ex: rotated_vect = R * vect.
%
% Outputs:
%   axis:       - Axis of rotation [x,y,z]
%   angle:      - Angle of rotation (given in radians)
% Methods from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/

%% Test for singularities:
% Symmetric matrices are associated with angle = 0 or pi
tolerance = 0.001;
if abs(R(1,2)-R(2,1))<tolerance && abs(R(1,3)-R(3,1))<tolerance && abs(R(3,2)-R(2,3))<tolerance
    %Matrix is considered symmetric
    if abs(R(1,2)+R(2,1))<tolerance && abs(R(1,3)+R(3,1))<tolerance && abs(R(3,2)+R(2,3))<tolerance
        % Only zeros outside diagonal
        if R(1,1) > 1-tolerance  && R(2,2) > 1-tolerance
            % 1 on the diagnonal
            angle = 0;
            % The axis doesn't matter
            axis = [1,0,0];
        else
            % one 1 and two -1 on the diagnonal
            angle = pi();
            if R(1,1) > 1-tolerance
                axis = [1,0,0];
            elseif R(2,2) > 1-tolerance
                axis = [0,1,0];
            else
                axis = [0,0,1];
            end
        end
    else
        angle = pi();
        % Find the highest diagonal term to minimize error
        dim = 1;
        dims2test = [2,3];
        for d = dims2test
            if R(dim,dim) < R(d,d)
                dim = d;
            end
        end
        
        axis = zeros(1,3);
        axis(dim) = sqrt((R(dim,dim)+1)/2);
        otherdims = setdiff(1:3, dim);
        for d = otherdims
            axis(d) = (R(dim,d)+R(d,dim))/(4*axis(dim));
        end
    end
else
    angle = acos((R(1,1)+R(2,2)+R(3,3)-1)/2);
    axis = zeros(1,3);
    divisor = sqrt((R(2,1)-R(1,2))^2 + (R(2,3)-R(3,2))^2 + (R(3,1)-R(1,3))^2);
    axis(1) = (R(3,2)-R(2,3))/divisor;
    axis(2) = (R(1,3)-R(3,1))/divisor;
    axis(3) = (R(2,1)-R(1,2))/divisor;
end
end