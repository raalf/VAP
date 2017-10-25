function [matVLST, matCENTER, matROTAX, matNEWWAKE, matNPNEWWAKE, vecDVEHVSPN, ...
    vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM] = ...
    fcnMOVEMULTIROTOR(vecRPM, matROTAX, matVLST, matCENTER, matNPVLST, ...
    vecDVETE, matDVE, vecAZNUM, valJ, valDIA, valALPHAR, vecDVEROTOR, vecDVEVLSTROTOR, vecRODIR)
% This function moves the rotor and created new wake elements in both
tempV = valJ*valDIA*vecRPM(1)/60;
vecJ = tempV./((vecRPM/60).*valDIA);

% Rotor tranlsation matrix
vecTRANSVLST = [-1*((vecJ(vecDVEVLSTROTOR).*valDIA)./vecAZNUM(vecDVEVLSTROTOR)).*cos(valALPHAR) zeros(size(vecDVEVLSTROTOR,1),1) ((vecJ(vecDVEVLSTROTOR).*valDIA)./vecAZNUM(vecDVEVLSTROTOR)).*sin(valALPHAR)];



tempVLST = matVLST - matROTAX(vecDVEVLSTROTOR,:);
tempCENTER = matCENTER - matROTAX(vecDVEROTOR,:);

% % Calculate momentum theory for wake displacement
% tempDELTIME = valAZNUM/(valRPM/60);
% tempDISPLACE = tempDELTIME*sqrt(-0.5.*(matUINF.^2) + sqrt((0.5.*(matUINF.^2)).^2+(mean(vecCTCONV)./(2*valDENSITY*(pi*(valDIA/2)^2)).^2)));

% Old trailing edge
matNEWWAKE(:,:,4) = matVLST(matDVE(vecDVETE>0,4),:);%+tempDISPLACE;
matNEWWAKE(:,:,3) = matVLST(matDVE(vecDVETE>0,3),:);%+tempDISPLACE;

% Old non-planer trailing edge
matNPNEWWAKE(:,:,4) = matNPVLST(matDVE(vecDVETE>0,4),:);%+tempDISPLACE;
matNPNEWWAKE(:,:,3) = matNPVLST(matDVE(vecDVETE>0,3),:);%+tempDISPLACE;

for i = 1:size(matVLST,1)
    valDELROT = vecRODIR(vecDVEVLSTROTOR(i)).*(2*pi)./vecAZNUM(vecDVEVLSTROTOR(i));
    % Rotor rotation matrix
    matROTATE = [cos(valDELROT) -sin(valDELROT) 0; sin(valDELROT) cos(valDELROT) 0; 0 0 1];

    % Move verticies and control points
    matVLST(i,:) = (matROTATE*(tempVLST(i,:))')' + matROTAX(vecDVEVLSTROTOR(i),:) + vecTRANSVLST(i,:);
end

for i = 1:size(matCENTER,1)
    valDELROT = vecRODIR(vecDVEROTOR(i))*(2*pi)/vecAZNUM(vecDVEROTOR(i));
    % Rotor rotation matrix
    matROTATE = [cos(valDELROT) -sin(valDELROT) 0; sin(valDELROT) cos(valDELROT) 0; 0 0 1];

    matCENTER(i,:) = (matROTATE*(tempCENTER(i,:))')' + matROTAX(vecDVEROTOR(i),:) + vecTRANSVLST(i,:);
end
% Define new trailing edge
matNEWWAKE(:,:,1) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,2) = matVLST(matDVE(vecDVETE>0,3),:);

% Define new non-planer trailing edge
matNPNEWWAKE(:,:,4) = matNPVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNPVLST(matDVE(vecDVETE>0,3),:);


% Move rotational axis (only translation)
tempTRANS = [-1*((vecJ.*valDIA)./vecAZNUM).*cos(valALPHAR) zeros(size(vecJ,1),1) ((vecJ.*valDIA)./vecAZNUM).*sin(valALPHAR)];
matROTAX = matROTAX + tempTRANS;

[ vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, ~, ~, ~ ] ...
    = fcnVLST2DVEPARAM(matDVE, matVLST);

end

