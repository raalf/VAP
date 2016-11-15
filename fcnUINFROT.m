function [matUINF] = fcnUINFROT(matCENTER, vecROTAX, valTIMESTEP, valRPM, valALPHAR, valAZNUM, valDIA, valJ)

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

% Find the current rotor angle (relative to initial position)
tempTHETA = valDELROT*valTIMESTEP;

% Convert rpm to rad per second
tempRADPS = valRPM*2*pi/60;

% Radial points of control points from rotational axis
tempCENTER = matCENTER - repmat(vecROTAX, length(matCENTER(:,1)),1);

% Radius Magnidute
tempRMAG = sqrt(tempCENTER(:,1).^2+tempCENTER(:,2).^2+tempCENTER(:,3).^2);

vecUROT = tempRADPS*tempRMAG*[cos(tempTHETA) sin(tempTHETA) 0];
vecUTRANS = (valJ*valDIA*valRPM/60)*[cos(valALPHAR) 0 sin(valALPHAR)];

matUINF = vecUROT + vecUTRANS;
end

