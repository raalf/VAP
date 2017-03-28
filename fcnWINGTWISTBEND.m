function [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt, vecLSM,...
    vecN, valSPAN, vecDVEHVSPN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE, valDELTIME, tempTIME, matDEF, matTWIST)
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

valDY = sum(2*vecDVEHVSPN,1)/length(vecDVEHVSPN);

valNSELE = sum(vecN,1)+1;

valSTRUCTDELTIME = valSDELTIME;

vecDEF = zeros(1,valNSELE+4);
vecTWIST = zeros(1,valNSELE+4);
vecSLOPE = zeros(1,valNSELE-1);

% Temporary cross sectional area calculation
C = -0.0333*vecSPANDIST + 0.76*ones(length(vecSPANDIST),1) ;
tk = 0.02 ;
Tk = 0.13 ;
vecSPANAREA = pi*tk*C*(1 + Tk);

valSTRUCTTIME = valTIMESTEP;
% valSTRUCTTIME = tempTIME + 2;

% vecJT = 0.0000001214.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST -  0.0000017210.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     + 0.0000051317.*vecSPANDIST.*vecSPANDIST - 0.0000073047.*vecSPANDIST + 0.0001181334;

% vecJT = 0.00045702.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     - 0.01320713.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     + 0.14939498.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST - 0.83266230.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     + 2.35858637.*vecSPANDIST.*vecSPANDIST - 3.18488527.*vecSPANDIST + 1.89798213;

% vecJT = 0.00004993.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     - 0.0015111.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     + 0.01788802.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST - 0.10454044.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
%     + 0.3130793.*vecSPANDIST.*vecSPANDIST - 0.45287673.*vecSPANDIST + 0.26571032;


%% Beam boundary conditions

matDEF(1:valSTRUCTTIME-1,:) = matDEFGLOB(1:valTIMESTEP-1,:);
% if tempTIME == 1
% matDEF(1:valSTRUCTTIME-1,:) = matDEFGLOB((valTIMESTEP-2):valTIMESTEP-1,:);
matTWIST(1:valSTRUCTTIME-1,:) = matTWISTGLOB(1:valTIMESTEP-1,:);
% matTWIST(1:valSTRUCTTIME-1,:) = matTWISTGLOB((valTIMESTEP-2):valTIMESTEP-1,:);
% end

vecDEF(3) = 0; % Zero deflection at root BC
vecTWIST(3) = 0; % Zero twist at root BC

% Assemble load matrix
matLOAD = [vecLIFTDIST' - vecLM.*9.81, vecMOMDIST - vecLM.*vecLSM.*9.81];

for yy = 4:(valNSELE+2)

    %% Geometric property assembly

    % Assemble mass matrix
%             matMASS = [vecLM(yy-2), -vecLM(yy-2).*vecLSM(yy-2); -vecLM(yy-2).*vecLSM(yy-2), vecLM(yy-2)*(vecLSM(yy-2)*vecLSM(yy-2) + vecJT(yy-2)/vecSPANAREA(yy-2))];
    matMASS = [vecLM(yy-2), -vecLM(yy-2).*vecLSM(yy-2); -vecLM(yy-2).*vecLSM(yy-2), vecJT(yy-2)];

    % Assemble stiffness matrices
    matK_1 = [matEIx(yy-2,3), 0; 0, 0];
    matK_2 = [matEIx(yy-2,2), 0; 0, -matGJt(yy-2,2)]; 
    matK_3 = [matEIx(yy-2,1), 0; 0, -matGJt(yy-2,1)];
    matB = [0 0; 0 3];

    %% Finite difference relations for partial derivatives

    % Finite difference relations for partial derivatives w.r.t
    % time
%     if tempTIME ~=1
        valUDOT = (matDEF(valSTRUCTTIME-1,yy) - matDEF(valSTRUCTTIME - 2, yy))./valSTRUCTDELTIME;
        valTDOT = (matTWIST(valSTRUCTTIME-1,yy) - matTWIST(valSTRUCTTIME - 2,yy))./valSTRUCTDELTIME;
%     else
%         valUDOT = (matDEF(valSTRUCTTIME-1,yy) - matDEF(valSTRUCTTIME - 2, yy))./valDELTIME;
%         valTDOT = (matTWIST(valSTRUCTTIME-1,yy) - matTWIST(valSTRUCTTIME - 2,yy))./valDELTIME;
%     end

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
    vecSLOPE(yy-3) = asin((vecDEF(yy)-vecDEF(yy-1))/(vecSPANDIST(yy-2)-vecSPANDIST(yy-3)));

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
% if tempTIME == floor(valDELTIME/valSDELTIME)
    matDEFGLOB(valTIMESTEP,:) = matDEF(end,:);
    matTWISTGLOB(valTIMESTEP,:) = matTWIST(end,:);

    matSLOPE(valTIMESTEP,:) = vecSLOPE';
% end

end