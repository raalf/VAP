function [matVLST, matCENTER, vecROTAX, matNEWWAKE, matNPNEWWAKE] = ...
    fcnMOVEROTOR(vecROTAX, valALPHAR, valAZNUM, valJ, matVLST, ...
    matCENTER, valDIA, matNPVLST, vecDVETE, matDVE)

% Finding the rotation angle for the timestep
valDELROT = (2*pi)/valAZNUM;

% Finding how much to translate the rotor based on alpha
vecL = [((valJ*valDIA)/valAZNUM)*cos(valALPHAR) 0 ((valJ*valDIA)/valAZNUM)*sin(valALPHAR)];


tempVLST = matVLST - repmat(vecROTAX, length(matVLST(:,1)),1);
tempCENTER = matCENTER - repmat(vecROTAX, length(matCENTER(:,1)),1);

% Defining rotation matrix
ROT = [cos(valDELROT) -sin(valDELROT) 0; sin(valDELROT) cos(valDELROT) 0; 0 0 1];

% Old trailing edge
matNEWWAKE(:,:,4) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,3) = matVLST(matDVE(vecDVETE>0,3),:);

% Old non-planer trailing edge
matNPNEWWAKE(:,:,4) = matNPVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNPVLST(matDVE(vecDVETE>0,3),:);

% Move verticies and control points
matVLST = (ROT*tempVLST')' + repmat(vecROTAX, length(matVLST(:,1)),1) + vecL;
matCENTER = (ROT*tempCENTER')' + repmat(vecROTAX, length(matCENTER(:,1)),1) + vecL;

% Define new trailing edge
matNEWWAKE(:,:,1) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,2) = matVLST(matDVE(vecDVETE>0,3),:);

% Define new non-planer trailing edge
matNPNEWWAKE(:,:,4) = matNPVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNPVLST(matDVE(vecDVETE>0,3),:);


% Move rotational axis (only translation)
vecROTAX = vecROTAX + vecL;

end

