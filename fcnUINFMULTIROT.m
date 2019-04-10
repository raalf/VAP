function [matUINF, matUINFTE, matTEPTS, vecTHETA] = fcnUINFMULTIROT(matCENTER, matROTAX, valTIMESTEP, vecRPM, valALPHAR, vecAZNUM, valDIA, valJ, valNUMB, vecDVEROTOR, vecDVEHVSPN, vecDVETE, matVLST, matDVE, vecRODIR)

% This function defines the direction and magnitude of the inflow velocity.
% This is calculated at both the control points and TE points.
%
% Note for TE matrix: (:,:,1) - left 80% halfspan of TE
%                     (:,:,2) - center of TE
%                     (:,:,3) - right 80% halfspan of TE
%
% OUTPUT:
%   matUINF - Velocity at each control point
%   matUINFTE - Velocity at each TE point of interest (for force
%   calculation)
%   matTEPTS - Trailing edge points of interest
%   vecTHETA - Current angle of each DVE
%   vecCPRADI - Radius from rotation axis to control point

tempV = valJ*valDIA*vecRPM(1)/60;
vecJ = tempV./((vecRPM/60).*valDIA);

% Finding the rotation angle for the timestep
valDELROT = (2.*pi)./vecAZNUM(vecDVEROTOR);

% Number of DVEs per blade
tempTOTDVE = size(matCENTER,1);
tempNUMDVE = tempTOTDVE/valNUMB;

% Find the current rotor angle (relative to initial position)
temp = 0:2*pi/(valNUMB):2*pi;
temp(valNUMB+1) = [];
vecTHETA = valDELROT.*valTIMESTEP.*vecRODIR(vecDVEROTOR)+reshape(repmat(temp,tempNUMDVE,1),tempTOTDVE,1);

% Convert rpm to rad per second
tempRADPS = vecRPM*2*pi/60;

% Radial points of control points from rotational axis
tempCENTER = matCENTER - matROTAX(vecDVEROTOR,:);

% Radius Magnitude
vecCPRADI = sqrt(tempCENTER(:,1).^2+tempCENTER(:,2).^2+tempCENTER(:,3).^2);

% Rotation matrix and traslation vector
matUROT =  repmat(tempRADPS(vecDVEROTOR),[1,3]).*[vecCPRADI.*cos(vecTHETA) vecCPRADI.*sin(vecTHETA) zeros(tempTOTDVE,1)];
vecUTRANS = (valDIA.*vecJ(vecDVEROTOR).*vecRPM(vecDVEROTOR)./60)*[-cos(valALPHAR) 0 sin(valALPHAR)];

if valJ == 0
    %vecUTRANS = ones(size(matCENTER,1)
end


% Velocity matrix for each control point
matUINF = matUROT - vecUTRANS;


%% Calculate trailing edge velocities
% Trailing edge points are calculated at %80 half span to the left and
% right of the TE as well as the center. Method is done as completed in
% fcnDVEINDDRAG of VAP.

% TE logical
idte = (vecDVETE == 3);

% Number of TE elements
numte = sum(idte);

% 80% of halfspan of TE elements only
eta8 = vecDVEHVSPN(idte).*0.8;

%TE vectors of TE elements only
s =( matVLST(matDVE(idte,3),:) -matVLST(matDVE(idte,4),:) )  ./ repmat((vecDVEHVSPN(idte).*2),1,3); %?? why is the S vector non-dim. to the span?

% Element TE edge midpoint of TE elements only
xte = (matVLST(matDVE(idte,3),:) + matVLST(matDVE(idte,4),:))/2;

matTEPTS = zeros(numte,3,3);

% Find 3 points along TE of all TE elements
% first layer is left side, second layer is middle, third layer is right
% side
matTEPTS(:,:,1) = (xte + s.*repmat(-eta8,1,3)); %left side
matTEPTS(:,:,2) = xte ; %middle
matTEPTS(:,:,3) = (xte + s.*repmat(eta8,1,3)); %right ride

% Radial points of trailing edge points from rotational axis
idxte = vecDVETE == 3;
tempTE = matTEPTS - repmat(matROTAX(vecDVEROTOR(idxte),:),1,1,3);

% Radius Magnitude
tempRMAGTE = sqrt(tempTE(:,1,:).^2+tempTE(:,2,:).^2+tempTE(:,3,:).^2);

% TE rotation matrix
matUROTTE = repmat(tempRADPS(vecDVEROTOR(idxte)),[1,3,3]).*[tempRMAGTE.*repmat(cos(vecTHETA(idte)),[1,1,3]) tempRMAGTE.*repmat(sin(vecTHETA(idte)),[1,1,3]) zeros(numte,1,3)];

% TE velocity matrix
matUINFTE = matUROTTE - repmat(vecUTRANS(idxte,:),[1,1,3]);

%quiver3(matTEPTS(:,1,2),matTEPTS(:,2,2),matTEPTS(:,3,2),matUINFTE(:,1,2),matUINFTE(:,2,2),matUINFTE(:,3,2))
%quiver3(matCENTER(:,1,1),matCENTER(:,2,1),matCENTER(:,3,1),matUINF(:,1,1),matUINF(:,2,1),matUINF(:,3,1))
end

