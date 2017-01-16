function [matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
    vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matVLST, matDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW)
%% Geometric Properties

% Find DVEs on LE and TE of wing
[ledves, ~, ~] = find(vecDVELE > 0);
[tedves, ~, ~] = find(vecDVETE > 0);

lepanels = vecDVEPANEL(ledves);

for i = 1:max(vecDVEWING)

    idxdve = ledves(vecDVEWING(ledves) == i);
    idxpanel = lepanels(vecDVEWING(ledves) == i);

    m = vecM(idxpanel);
    if any(m - m(1))
        disp('Problem with wing chordwise elements.');
        break
    end
    m = m(1);

    tempm = repmat(vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel),1);

    rows = repmat(idxdve,1,m) + tempm; % DVEs along each chord station

end

tempDVEEDGECRD = abs(matVLST(matDVE(ledves,1),:) - matVLST(matDVE(tedves,4),:));
tempDVEEDGECRD = [tempDVEEDGECRD; abs(matVLST(matDVE(ledves(end),2),:) - matVLST(matDVE(tedves(end),3),:))];

tempDVEEDGECRD = sqrt(sum(tempDVEEDGECRD.^2,2)); % Chord length at each DVE edge

% Use transformation matrix to determine X,Y,Z coordinates of aerodynamic
% center based on DVE edge chord
matDVEEDGECRD = [tempDVEEDGECRD, zeros(length(tempDVEEDGECRD),2)]; % Chord vector at each DVE edge in local DVE frame (vector pointing from LE to TE)

matQTRCRD = fcnSTARGLOB(0.25*matDVEEDGECRD,[vecDVEROLL(rows(:,1));vecDVEROLL(rows(end,1))],[vecDVEPITCH(rows(:,1));vecDVEPITCH(rows(end,1))],[vecDVEYAW(rows(:,1));vecDVEYAW(rows(end,1))]);

tempLE = [matVLST(matDVE(ledves,1),:); matVLST(matDVE(ledves(end),2),:)]; % DVE LE coordinates

matAEROCNTR = tempLE + matQTRCRD; % X, Y, Z location of aerodynamic center

% Determine spanwise location (y coordinate) of DVEs
tempSPANDIST = 2*vecDVEHVSPN(ledves);
matSPANDIST = repmat(tempSPANDIST,1,size(tempSPANDIST));

tempSPAN = triu(matSPANDIST);  % Upper triangular matrix of DVE spans

vecSPANDIST = sum(tempSPAN(:,1:size(matSPANDIST,2)))';
vecSPANDIST = [0; vecSPANDIST]; % Adding location of wing root


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

% Determining X,Y,Z location of elastic axis (shear center)
tempEA = [vecEA.*tempDVEEDGECRD, zeros(length(vecEA),2)] + matDVEEDGECRD*0.25; % Distance to EA from LE in local coordinates

matSC = fcnSTARGLOB(tempEA,[vecDVEROLL(rows(:,1));vecDVEROLL(rows(end,1))],[vecDVEPITCH(rows(:,1));vecDVEPITCH(rows(end,1))],[vecDVEYAW(rows(:,1));vecDVEYAW(rows(end,1))]); % Transform to global coordinates

matSC = tempLE + matSC; % Add LE coordinates to have absolute location

matLSAC = matAEROCNTR - matSC;

vecLSAC = sqrt(sum(matLSAC.^2,2));

vecLSM = 0.1.*vecLSAC;

% Determine distance of each vertex to SC (to be used for twist velocity
% later)
temp_leftV = [matDVE(rows,1),matDVE(rows,4)];
temp_leftV = reshape(temp_leftV,sum(vecN,1),[]);

[move_row,~] = find(temp_leftV); % Vector corresponding to which shear center coordinate to use

tempSCLST = zeros(size(matVLST,1),3);

tempSCLST(temp_leftV,:) = matSC(move_row,:);

temp_rightV = [matDVE(rows,2), matDVE(rows,3)];
temp_rightV = reshape(temp_rightV,sum(vecN,1),[]);

[move_row,~] = find(temp_rightV); % Vector corresponding to which shear center coordinate to use
tempSCLST(temp_rightV,:) = matSC(move_row,:);

matSCLST = tempSCLST;

matSCLST = matSCLST - matVLST; % Matrix of vectors between shear center and vertex

end


