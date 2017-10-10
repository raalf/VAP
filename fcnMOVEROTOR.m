function [matVLST, matCENTER, vecROTAX, matNEWWAKE, matNPNEWWAKE, vecDVEHVSPN, ...
    vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM] = ...
    fcnMOVEROTOR(valDENSITY, valRPM, vecROTAX, vecCTCONV, matUINF, matVLST, matCENTER, matNPVLST, ...
    vecDVETE, matDVE, valAZNUM, valJ, valDIA, valALPHAR)
% This function moves the rotor and created new wake elements in both
% tranlational and rotational directions. 

valDELROT = (2*pi)/valAZNUM;

% Rotor rotation matrix
matROTATE = [cos(valDELROT) -sin(valDELROT) 0; sin(valDELROT) cos(valDELROT) 0; 0 0 1];

% Rotor tranlsation matrix
vecTRANS = [-1*((valJ*valDIA)/valAZNUM)*cos(valALPHAR) 0 ((valJ*valDIA)/valAZNUM)*sin(valALPHAR)];

tempVLST = matVLST - repmat(vecROTAX, length(matVLST(:,1)),1);
tempCENTER = matCENTER - repmat(vecROTAX, length(matCENTER(:,1)),1);

% % Calculate momentum theory for wake displacement
% tempDELTIME = valAZNUM/(valRPM/60);
% tempDISPLACE = tempDELTIME*sqrt(-0.5.*(matUINF.^2) + sqrt((0.5.*(matUINF.^2)).^2+(mean(vecCTCONV)./(2*valDENSITY*(pi*(valDIA/2)^2)).^2)));

% Old trailing edge
matNEWWAKE(:,:,4) = matVLST(matDVE(vecDVETE>0,4),:);%+tempDISPLACE;
matNEWWAKE(:,:,3) = matVLST(matDVE(vecDVETE>0,3),:);%+tempDISPLACE;

% Old non-planer trailing edge
matNPNEWWAKE(:,:,4) = matNPVLST(matDVE(vecDVETE>0,4),:);%+tempDISPLACE;
matNPNEWWAKE(:,:,3) = matNPVLST(matDVE(vecDVETE>0,3),:);%+tempDISPLACE;

% Move verticies and control points
matVLST = (matROTATE*tempVLST')' + repmat(vecROTAX, length(matVLST(:,1)),1) + vecTRANS;
matCENTER = (matROTATE*tempCENTER')' + repmat(vecROTAX, length(matCENTER(:,1)),1) + vecTRANS;

% Define new trailing edge
matNEWWAKE(:,:,1) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,2) = matVLST(matDVE(vecDVETE>0,3),:);

% Define new non-planer trailing edge
matNPNEWWAKE(:,:,4) = matNPVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNPVLST(matDVE(vecDVETE>0,3),:);


% Move rotational axis (only translation)
vecROTAX = vecROTAX + vecTRANS;

[ vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, ~, ~, ~ ] ...
    = fcnVLST2DVEPARAM(matDVE, matVLST);

end

