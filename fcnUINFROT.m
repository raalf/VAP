function [matUINF] = fcnUINFROT(matCENTER, vecROTAX, valTIMESTEP, valRPM, valALPHAR, valAZNUM, valDIA, valJ, valNUMB)

% This function defines the direction of the inflow velocity

% INPUT:
%   valALPHAR - Angle of attack (radians)
%   valAXNUM - Number of azmith locations to find timestep distance
%   valDIA - Rotor diameter
%   valJ - Advance ratio
%
% OUTPUT:
%   matROTATE - Rotation matrix
%   vecTRANS - Translation matrix   


% Finding the rotation angle for the timestep
valDELROT = (2*pi)/valAZNUM;

% Number of DVEs per blade
tempTOTDVE = size(matCENTER);
tempNUMDVE = tempTOTDVE(1)/valNUMB;

% Find the current rotor angle (relative to initial position)
temp = 0:2*pi/(valNUMB):2*pi;
temp(valNUMB+1) = [];
vecTHETA = valDELROT*valTIMESTEP+reshape(repmat(temp,tempNUMDVE,1),tempTOTDVE(1),1);


% Convert rpm to rad per second
tempRADPS = valRPM*2*pi/60;

% Radial points of control points from rotational axis
tempCENTER = matCENTER - repmat(vecROTAX, tempTOTDVE(1),1);

% Radius Magnitude
tempRMAG = sqrt(tempCENTER(:,1).^2+tempCENTER(:,2).^2+tempCENTER(:,3).^2);

matUROT = tempRADPS*[tempRMAG.*cos(vecTHETA) tempRMAG.*sin(vecTHETA) zeros(tempTOTDVE(1),1)];
vecUTRANS = (-1*valJ*valDIA*valRPM/60)*[cos(valALPHAR) 0 sin(valALPHAR)];

matUINF = matUROT - vecUTRANS;


%quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),matUINF(:,1),matUINF(:,2),matUINF(:,3))
end

