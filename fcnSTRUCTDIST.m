function [matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST, matSC, vecMAC] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
    vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW)
%% Geometric Properties

% Find DVEs on LE and TE of wing
[ledves, ~, ~] = find(vecDVELE > 0);
[tedves, ~, ~] = find(vecDVETE > 0);

[matROWS] = fcnDVEROW(ledves, vecDVEPANEL, vecDVEWING, vecM, vecN);

tempDVEEDGECRD = abs(matNPVLST(matNPDVE(ledves,1),:) - matNPVLST(matNPDVE(tedves,4),:));
tempDVEEDGECRD = [tempDVEEDGECRD; abs(matNPVLST(matNPDVE(ledves(end),2),:) - matNPVLST(matNPDVE(tedves(end),3),:))];

tempDVEEDGECRD = sqrt(sum(tempDVEEDGECRD.^2,2)); % Chord length at each DVE edge

% Use transformation matrix to determine X,Y,Z coordinates of aerodynamic
% center based on DVE edge chord
matDVEEDGECRD = [tempDVEEDGECRD, zeros(length(tempDVEEDGECRD),2)]; % Chord vector at each DVE edge in local DVE frame (vector pointing from LE to TE)

matQTRCRD = fcnSTARGLOB(0.25*matDVEEDGECRD,[vecDVEROLL(matROWS(:,1));vecDVEROLL(matROWS(end,1))],[vecDVEPITCH(matROWS(:,1));vecDVEPITCH(matROWS(end,1))],[vecDVEYAW(matROWS(:,1));vecDVEYAW(matROWS(end,1))]);

tempLE = [matNPVLST(matNPDVE(ledves,1),:); matNPVLST(matNPDVE(ledves(end),2),:)]; % DVE LE coordinates

matAEROCNTR = tempLE + matQTRCRD; % X, Y, Z location of aerodynamic center

% Determine spanwise location (y coordinate) of DVEs
tempSPANDIST = 2*vecDVEHVSPN(ledves);
matSPANDIST = repmat(tempSPANDIST,1,size(tempSPANDIST));

tempSPAN = triu(matSPANDIST);  % Upper triangular matrix of DVE spans

vecSPANDIST = sum(tempSPAN(:,1:size(matSPANDIST,2)))';
vecSPANDIST = [0; vecSPANDIST]; % Adding location of wing root

% Calculating mean aerodynamic chord at each spanwise station to use for
% pitching moment calculations later on

idx1 = 1:(length(tempDVEEDGECRD)-1); % Index of root chord elements
idx2 = 2:length(tempDVEEDGECRD); % Index of tip chord elements

% Calculate taper ratio of each DVE
tempDVEEDGECRD = repmat(tempDVEEDGECRD,1,2);
taper_ratio = (tempDVEEDGECRD(idx2)./tempDVEEDGECRD(idx1))';

vecMAC = (2/3)*tempDVEEDGECRD(idx1,1).*(1 + taper_ratio + taper_ratio.^2)./(1 + taper_ratio); % Vector of mean aerodynamic chord at each spanwise station


%% Structural Properties

% Spanwise bending stiffness distribution. Cols 2 and 3 are the first and
% second derivatives
matEIx(:,1) = (vecEIxCOEFF(1).*vecSPANDIST.^2 + vecEIxCOEFF(2).*vecSPANDIST + vecEIxCOEFF(3))';
matEIx(:,2) = (2*vecEIxCOEFF(1).*vecSPANDIST + vecEIxCOEFF(2))';
matEIx(:,3) = (repmat(2*vecEIxCOEFF(1),1,size(vecSPANDIST)))';

% Spanwise torsional stiffness distribution. Col 2 is the first derivative
matGJt(:,1) = (vecGJtCOEFF(1).*vecSPANDIST.^2 + vecGJtCOEFF(2).*vecSPANDIST + vecGJtCOEFF(3))';
matGJt(:,2) = (2*vecGJtCOEFF(1).*vecSPANDIST + vecGJtCOEFF(2))';

vecEA = vecEACOEFF(1).*vecSPANDIST.^2 + vecEACOEFF(2).*vecSPANDIST + vecEACOEFF(3);
vecCG = vecCGCOEFF(1).*vecSPANDIST.^2 + vecCGCOEFF(2).*vecSPANDIST + vecCGCOEFF(3);
vecJT = vecJTCOEFF(1).*vecSPANDIST.^2 + vecJTCOEFF(2).*vecSPANDIST + vecJTCOEFF(3);
vecLM = vecLMCOEFF(1).*vecSPANDIST.^2 + vecLMCOEFF(2).*vecSPANDIST + vecLMCOEFF(3);

% vecLM(end) = 0; % Add zero mass to wing tip

% Determining X,Y,Z location of elastic axis (shear center) and center of
% mass (CG)
tempEA = [vecEA, zeros(length(vecEA),2)]; % Distance to EA from LE in local coordinates

tempCG = [vecCG, zeros(length(vecCG),2)]; % Distance to CG from LE in local coordinates

matCG = fcnSTARGLOB(tempCG,[vecDVEROLL(matROWS(:,1));vecDVEROLL(matROWS(end,1))],[vecDVEPITCH(matROWS(:,1));vecDVEPITCH(matROWS(end,1))],[vecDVEYAW(matROWS(:,1));vecDVEYAW(matROWS(end,1))]); % Transform to global coordinates

matSC = fcnSTARGLOB(tempEA,[vecDVEROLL(matROWS(:,1));vecDVEROLL(matROWS(end,1))],[vecDVEPITCH(matROWS(:,1));vecDVEPITCH(matROWS(end,1))],[vecDVEYAW(matROWS(:,1));vecDVEYAW(matROWS(end,1))]); % Transform to global coordinates

matSC = tempLE + matSC; % Add LE coordinates to have absolute location

matCG = tempLE + matCG;

matLSM = matSC - matCG;

matLSAC = matAEROCNTR - matSC;

vecLSAC = sqrt(sum(matLSAC.^2,2));

vecLSM = -1.*sign(matLSM(:,1)).*sqrt(sum(matLSM.^2,2)); % If +ve --> CG is ahead of SC; If -ve --> CG is behind SC

% Determine distance of each vertex to SC (to be used for twist velocity
% later)
temp_leftV = [matNPDVE(matROWS,1),matNPDVE(matROWS,4)];
temp_leftV = reshape(temp_leftV,sum(vecN,1),[]);

[move_row,~] = find(temp_leftV); % Vector corresponding to which shear center coordinate to use

tempSCLST = zeros(size(matNPVLST,1),3);

tempSCLST(temp_leftV,:) = matSC(move_row,:);

temp_rightV = [matNPDVE(matROWS,2), matNPDVE(matROWS,3)];
temp_rightV = reshape(temp_rightV,sum(vecN,1),[]);

[move_row,~] = find(temp_rightV); % Vector corresponding to which shear center coordinate to use
tempSCLST(temp_rightV,:) = matSC(move_row+1,:);

matSCLST = tempSCLST;

matSCLST = matSCLST - matNPVLST; % Matrix of vectors between shear center and vertex

% figure(4)
% clf
% patch('Faces',matNPDVE,'Vertices',matNPVLST,'FaceColor','r')
% hold on
% plot3(matAEROCNTR(:,1), matAEROCNTR(:,2), matAEROCNTR(:,3),'-ok')
% plot3(matSC(:,1), matSC(:,2), matSC(:,3),'-ob')

end


