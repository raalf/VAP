function [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND_STAGGER(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt, vecLSM,...
    valSPAN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE, matDEF, matTWIST, valNSELE, tempTIME)
% This function computes the spanwise deflection and twist using an
% explicit finite difference method given a loading and structural
% distribution.
%
% INPUTS:
%
% valTIMESTEP - Current timestep number
%
% vecLIFTDIST - 1 x n vector with the lift values at each node, where n is
% the number of spanwise stations
%
% vecMOMDIST - 1 x n vector with aerodynamic moment at each node, where n
% is the number of spanwise stations
%
% vecSPANAREA - 1 x n vector containing the structural cross sectional area
% at each spanwise location, where n is the number of spanwise stations
%
% matEIx - n x 3 matrix containing the bending stiffness at each spanwise
% station, where n is the number of spanwise stations. The first column
% represents EIx, the second column EIx', and the third column EIx''
%
% matGJt - n x 2 matrix containing the torsional stiffness distribution at
% each spanwise node, where n is the number of spanwise stations. The first
% row is GJt and the secon GJt'
%
% vecTORSIONRIGIDITY - 1 x n vector containing the torsional rigidity (GJ)
% at each spanwise station, where n is the number of spanwise stations.
%
% vecMASS2SHEAR - 1 x n vector containing the distances (in m) between the
% center of mass and shear center at each spanwise station, where n is the
% number of spanwise stations
%
% valYMODULUS - Young's modulus of the wing structure in Pa
%
% valNSELE - Number of spanwise elements
%

valDY = 0.5*valSPAN/valNSELE;

temp_y = (0:valDY:0.5*valSPAN)';

vecSPANDIST(end) = temp_y(end);

% valNSELE = sum(vecN,1)+1;

valSTRUCTDELTIME = valSDELTIME;

% vecJT = 0.00037078.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     - 0.01102270.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     + 0.12838255.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST - 0.73708913.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     + 2.15067037.*vecSPANDIST.*vecSPANDIST - 2.99312818.*vecSPANDIST + 1.84576176;

%% Interpolate values at structural grid points

[matEIx_interp(:,1)] = interp1(vecSPANDIST,matEIx(:,1)',temp_y);
[matEIx_interp(:,2)] = interp1(vecSPANDIST,matEIx(:,2)',temp_y);
[matEIx_interp(:,3)] = interp1(vecSPANDIST,matEIx(:,3)',temp_y);
[matGJt_interp(:,1)] = interp1(vecSPANDIST,matGJt(:,1)',temp_y);
[matGJt_interp(:,2)] = interp1(vecSPANDIST,matGJt(:,2)',temp_y);
[vecLM] = interp1(vecSPANDIST,vecLM',temp_y);
[vecJT] = interp1(vecSPANDIST,vecJT',temp_y);
[vecLSM] = interp1(vecSPANDIST,vecLSM',temp_y);
[vecLIFTDIST] = interp1(vecSPANDIST,vecLIFTDIST,temp_y);
[vecMOMDIST] = interp1(vecSPANDIST,vecMOMDIST',temp_y);

% vecLSM = -vecLSM;
matEIx = matEIx_interp;
matGJt = matGJt_interp;

vecDEF = zeros(1,valNSELE+4);
vecTWIST = zeros(1,valNSELE+4);
vecSLOPE = zeros(1,valNSELE-1);

% valSTRUCTTIME = valTIMESTEP;
valSTRUCTTIME = tempTIME + 2;

%% Beam boundary conditions

% matDEF(1:valSTRUCTTIME-1,:) = matDEFGLOB(1:valTIMESTEP-1,:);
% matDEF(1:valSTRUCTTIME-1,:) = matDEF(1:valTIMESTEP-1,:);
% if tempTIME == 1
% matDEF(1:valSTRUCTTIME-1,:) = matDEFGLOB((valTIMESTEP-2):valTIMESTEP-1,:);
% matTWIST(1:valSTRUCTTIME-1,:) = matTWISTGLOB(1:valTIMESTEP-1,:);
% matTWIST(1:valSTRUCTTIME-1,:) = matTWIST(1:valTIMESTEP-1,:);
% matTWIST(1:valSTRUCTTIME-1,:) = matTWISTGLOB((valTIMESTEP-2):valTIMESTEP-1,:);
% end

if tempTIME == 1
    matDEF(1:valSTRUCTTIME-1,:) = matDEF((end-1):end,:);
    matTWIST(1:valSTRUCTTIME-1,:) = matTWIST((end-1):end,:);
end

vecDEF(3) = 0; % Zero deflection at root BC
vecTWIST(3) = 0; % Zero twist at root BC

% Assemble load matrix
matLOAD = [vecLIFTDIST - vecLM.*9.81, vecMOMDIST + vecLSM.*vecLM.*9.81];
% matLOAD = [vecLIFTDIST', vecMOMDIST'];
% matLOAD(end,:) = [0,0]; 

for yy = 4:(valNSELE+2)

    %% Geometric property assembly

    % Assemble mass matrix
    matMASS = [vecLM(yy-2), -vecLM(yy-2).*vecLSM(yy-2); -vecLM(yy-2).*vecLSM(yy-2), vecJT(yy-2)];

    % Assemble stiffness matrices
    matK_1 = [matEIx(yy-2,3), 0; 0, 0];
    matK_2 = [matEIx(yy-2,2), 0; 0, -matGJt(yy-2,2)];
    matK_3 = [matEIx(yy-2,1), 0; 0, -matGJt(yy-2,1)];
    matB = [100 0; 0 200];

    %% Finite difference relations for partial derivatives

    % Finite difference relations for partial derivatives w.r.t
    % time
    valUDOT = (matDEF(valSTRUCTTIME-1,yy) - matDEF(valSTRUCTTIME - 2, yy))./valSTRUCTDELTIME;
    valTDOT = (matTWIST(valSTRUCTTIME-1,yy) - matTWIST(valSTRUCTTIME - 2,yy))./valSTRUCTDELTIME;

    % Finite difference relations for partial derivative of deflection w.r.t Y
    valU_yy = (matDEF(valSTRUCTTIME-1,yy+1) - 2*matDEF(valSTRUCTTIME-1,yy) + matDEF(valSTRUCTTIME-1,yy-1))/(valDY)^2;
    valU_yyy = (matDEF(valSTRUCTTIME-1,yy+2) - 3*matDEF(valSTRUCTTIME-1,yy+1) + 3*matDEF(valSTRUCTTIME-1,yy)- ...
        matDEF(valSTRUCTTIME-1,yy-1))/(valDY)^3;
    valU_yyyy = (matDEF(valSTRUCTTIME-1,yy+2) - 4*matDEF(valSTRUCTTIME-1,yy+1) + 6*matDEF(valSTRUCTTIME-1,yy) - ...
        4*matDEF(valSTRUCTTIME-1,yy-1) + matDEF(valSTRUCTTIME-1,yy-2))/(valDY)^4;

    % Finite difference relations for partial derivative of twist w.r.t Y
    valTHETA_y = (matTWIST(valSTRUCTTIME-1,yy+1) - matTWIST(valSTRUCTTIME-1,yy-1))/(2*valDY);
    valTHETA_yy = (matTWIST(valSTRUCTTIME-1,yy+1) - 2*matTWIST(valSTRUCTTIME-1,yy) + matTWIST(valSTRUCTTIME-1,yy-1))/(valDY^2);

    %% Solve matrix equation

    % Temp variable with the wing deflection and twist stored as a matrix. The
    % first row is the deflection, w/ each column as a spanwise station. The
    % second row is the twist, w/ each column as a spanwise station.

    tempTWISTBEND = 2.*[matDEF(valSTRUCTTIME-1,yy); matTWIST(valSTRUCTTIME-1,yy)] - [matDEF(valSTRUCTTIME-2,yy); matTWIST(valSTRUCTTIME-2,yy)] ...
        + (valSTRUCTDELTIME^2).*inv(matMASS)*([matLOAD(yy-2,1); matLOAD(yy-2,2)] - matK_1*[valU_yy; 0] - matK_2*[valU_yyy; valTHETA_y] -...
        matK_3*[valU_yyyy; valTHETA_yy] - matB*[valUDOT; valTDOT]);

    % Output result of deflection and twist to separate vectors
    vecDEF(yy) = tempTWISTBEND(1,:);
    vecTWIST(yy) = tempTWISTBEND(2,:);

    % Calculate angle between DVE and horizontal based on
    % deflection
    vecSLOPE(yy-3) = asin((vecDEF(yy)-vecDEF(yy-1))/(temp_y(yy-2)-temp_y(yy-3)));

end

vecDEF(valNSELE+3) = 2*vecDEF(valNSELE+2)...
    -vecDEF(valNSELE+1); % BC for deflection one element beyond wing (positive span direction)

vecDEF(valNSELE+4) = 3*vecDEF(valNSELE+2)...
    -2*vecDEF(valNSELE+1); % BC for deflection two elements beyond wing (positive span direction)

vecDEF(2) = vecDEF(4); % BC for deflection one element beyond root (negative span direction)

vecTWIST(valNSELE+3) = vecTWIST(valNSELE+1); % BC for twist one element beyond wing tip (positive span direction)

matDEF(valSTRUCTTIME,:) = vecDEF;
matTWIST(valSTRUCTTIME,:) = vecTWIST;

% Spanwise deflection and twist wrt structural timestep
vecDEF = matDEF(end,:);
vecTWIST = matTWIST(end,:);

vecSLOPE = [0; vecSLOPE'];

% Spanwise deflection and twist wrt to global timestep
matDEFGLOB(valTIMESTEP,:) = interp1(temp_y,matDEF(valSTRUCTTIME,3:end-1),vecSPANDIST);
matTWISTGLOB(valTIMESTEP,:) = interp1(temp_y,matTWIST(valSTRUCTTIME,3:end-1),vecSPANDIST);

matSLOPE(valTIMESTEP,:) = interp1(temp_y(1:end-1),vecSLOPE,vecSPANDIST(1:end-1));


end