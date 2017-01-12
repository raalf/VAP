function [matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, vecSPANDIST] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
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

vecLSAC = (0.0062*vecSPANDIST.*vecSPANDIST.*vecSPANDIST - 0.0533*vecSPANDIST.*vecSPANDIST + 0.1403*vecSPANDIST + 0.7029)*0.01;

test = 0.25*matDVEEDGECRD + [vecLSAC, zeros(length(vecLSAC),2)];

tempLSM = 0.1*test;

tempSC = fcnSTARGLOB(test,[vecDVEROLL(rows(:,1));vecDVEROLL(rows(end,1))],[vecDVEPITCH(rows(:,1));vecDVEPITCH(rows(end,1))],[vecDVEYAW(rows(:,1));vecDVEYAW(rows(end,1))]);

tempLSM2 = fcnSTARGLOB(tempLSM,[vecDVEROLL(rows(:,1));vecDVEROLL(rows(end,1))],[vecDVEPITCH(rows(:,1));vecDVEPITCH(rows(end,1))],[vecDVEYAW(rows(:,1));vecDVEYAW(rows(end,1))]);

matSC = tempLE + tempSC;

matLSM = tempLE + tempLSM2;


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

test = fcnSTARGLOB([test, zeros(length(test),2)],[vecDVEROLL(rows(:,1));vecDVEROLL(rows(end,1))],[vecDVEPITCH(rows(:,1));vecDVEPITCH(rows(end,1))],[vecDVEYAW(rows(:,1));vecDVEYAW(rows(end,1))]);

testSC = test + tempLE;

matSC2 = [vecEA, matAEROCNTR(:,2), matAEROCNTR(:,3)];

vecLSM = 0;

end


