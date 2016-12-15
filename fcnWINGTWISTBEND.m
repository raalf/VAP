function [vecDEFLECTION, vecTWIST] = fcnWINGTWISTBEND(vecLIFTDIST, vecMOMDIST, matDEFLECTION, matTWIST, vecSPANAREA, matEIx, vecLM, matGJt, vecLSM, valYMODULUS, vecN, valSPAN, vecDEFLECTION, vecTWIST, vecBEAM, valINITCOND, valDY)
%% Function Inputs
%
% valDELTIME - Timestep size (s)
%
% valTIMESTEP - Current timestep number
%
% vecLIFT - 1 x n vector with the lift values at each node, where n is the
% number of spanwise stations
%
% vecMOM - 1 x n vector with aerodynamic moment at each node, where n is
% the number of spanwise stations
%
% matDEFLECTION - valTIMESTEP x n matrix containing the deflection at each
% spanwise station. Rows are solutions at each previous timestep(s),
% columns are spanwise locations.
%
% matTWIST - valTIMESTEP x n matrix containing the twist at each spanwise
% station. Rows are solutions at previous timestep(s), columns are spanwise
% locations.
%
% vecSPANAREA - 1 x n vector containing the structural cross sectional area
% at each spanwise location, where n is the number of spanwise stations
%
% matSPANINERTIA - 3 x n matrix containing the area moment of inertia at
% each spanwise station, where n is the number of spanwise stations. The
% first row represents Ixx, the second row Ixx', and the third row Ixx''
%
% vecLINMASS - 1 x (n-1) vector containing the linear mass in kg/m of each
% element, where n is the number of spanwise stations.
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
valNSELE = size(vecN,1);

valTIMESTEP = 0;

valDELTIME = sqrt(150000*valDY^4/valYMODULUS);
%% Beam boundary conditions
   
if valTIMESTEP == 1

    vecDEFLECTION(3) = 0; % Zero deflection at root BC
    vecDEFLECTION(2) = vecDEFLECTION(4); % BC for deflection one element beyond root (negative span direction)
    vecDEFLECTION(1) = vecDEFLECTION(5); % BC for deflection two elements beyond root (negative span direction)
    
    vecDEFLECTION(3:valNSELE+2) = valINITCOND*vecBEAM.*vecBEAM;
    % vecDEFLECTION(3:valNSELE+2) = 0;
    
    vecDEFLECTION(valNSELE+3) = 2*vecDEFLECTION(valNSELE+2)...
        -vecDEFLECTION(valNSELE+1); % BC for deflection one element beyond wing (positive span direction)
    vecDEFLECTION(valNSELE+4) = 3*vecDEFLECTION(valNSELE+2)...
        -2*vecDEFLECTION(valNSELE+1); % BC for deflection two elements beyond wing (positive span direction)   
    
    vecTWIST(3) = 0;
    vecTWIST(2) = vecTWIST(4);
    vecTWIST(1) = vecTWIST(5);
    
    vecTWIST(3:valNSELE+2) = valINITCOND*(vecBEAM./valSPAN);
    
    vecTWIST(valNSELE+3) = vecTWIST(valNSELE+1); % BC for twist one element beyond wing (positive span direction)

    matDEFLECTION(valTIMESTEP,:) = vecDEFLECTION;
    matTWIST(valTIMESTEP,:) = vecTWIST;

elseif valTIMESTEP == 2

    vecDEFLECTION(3) = 0; % Zero deflection at root BC
    vecDEFLECTION(2) = vecDEFLECTION(4); % BC for deflection one element beyond root (negative span direction)
    vecDEFLECTION(1) = vecDEFLECTION(5); % BC for deflection two elements beyond root (negative span direction)
    
    vecDEFLECTION(3:valNSELE+2) = valINITCOND*vecBEAM.*vecBEAM;
    % vecDEFLECTION(3:valNSELE+2) = 0;
    
    vecDEFLECTION(valNSELE+3) = 2*vecDEFLECTION(valNSELE+2)...
        -vecDEFLECTION(valNSELE+1); % BC for deflection one element beyond wing (positive span direction)
    vecDEFLECTION(valNSELE+4) = 3*vecDEFLECTION(valNSELE+2)...
        -2*vecDEFLECTION(valNSELE+1); % BC for deflection two elements beyond wing (positive span direction)

    vecTWIST(3) = 0;
    vecTWIST(2) = vecTWIST(4);
    vecTWIST(1) = vecTWIST(5);
    
    vecTWIST(3:valNSELE+2) = valINITCOND*(vecBEAM./valSPAN);
    
    vecTWIST(valNSELE+3) = vecTWIST(valNSELE+1); % BC for twist one element beyond wing (positive span direction)

    matDEFLECTION(valTIMESTEP,:) = vecDEFLECTION;
    matTWIST(valTIMESTEP,:) = vecTWIST;

else

    for yy = 4:(valNSELE+2)
        %% Geometric property assembly

        % Assemble mass matrix
        matMASS = [vecLM(yy-2), -vecLM(yy-2).*vecLSM(yy-2); -vecLM(yy-2).*vecLSM(yy-2), vecLM(yy-2).*(vecLSM(yy-2)+(matGJt(yy-2)./(vecSPANAREA(yy-2))))];

        % Assemble load matrix
        matLOAD = [vecLIFTDIST(yy-2); vecMOMDIST(yy-2)];

        % Assemble stiffness matrices
        matK_1 = [valYMODULUS.*matEIx(3,yy-2), 0; 0, 0]; % Need to figure out second derivative of I (I'')... standy by
        matK_2 = [valYMODULUS.*matEIx(2,yy-2), 0; 0, 0]; % Need to figure out derivative of I (I')
        matK_3 = [valYMODULUS.*matEIx(1,yy-2), 0; 0, -matGJt(yy-2)];

        %% Finite difference relations for partial derivatives

        % Finite difference relations for partial derivative of deflection w.r.t Y
        valU_yy = (matDEFLECTION(valTIMESTEP-1,yy+1) - 2*matDEFLECTION(valTIMESTEP-1,yy) + matDEFLECTION(valTIMESTEP-1,yy-1))/(valDY)^2;
        valU_yyy = (matDEFLECTION(valTIMESTEP-1,yy+2) - 3*matDEFLECTION(valTIMESTEP-1,yy+1) + 3*matDEFLECTION(valTIMESTEP-1,yy)- ...
            matDEFLECTION(valTIMESTEP-1,yy-1))/(valDY)^3;
        valU_yyyy = (matDEFLECTION(valTIMESTEP-1,yy+2) - 4*matDEFLECTION(valTIMESTEP-1,yy+1) + 6*matDEFLECTION(valTIMESTEP-1,yy) - ...
            4*matDEFLECTION(valTIMESTEP-1,yy-1) + matDEFLECTION(valTIMESTEP-1,yy-2))/(valDY)^4;

        % Finite difference relations for partial derivative of twist w.r.t Y
        valTHETA_y = (matTWIST(valTIMESTEP-1,yy+1) - matTWIST(valTIMESTEP-1,yy-1))/2*valDY;
        valTHETA_yy = (matTWIST(valTIMESTEP-1,yy+1) - 2*matTWIST(valTIMESTEP-1,yy) + matTWIST(valTIMESTEP-1,yy-1))/(valDY^2);

        %% Solve matrix equation

        % Temp variable with the wing deflection and twist stored as a matrix. The
        % first row is the deflection, w/ each column as a spanwise station. The
        % second row is the twist, w/ each column as a spanwise station.

        vecDEFLECTION(3) = 0; % Zero deflection at root BC
        vecTWIST(3) = 0;

        tempTWISTBEND = 2.*[matDEFLECTION(valTIMESTEP-1,yy); matTWIST(valTIMESTEP-1,yy)] - [matDEFLECTION(valTIMESTEP-2,yy); matTWIST(valTIMESTEP-2,yy)] ...
            + (valDELTIME^2).*inv(matMASS)*(matLOAD - matK_1*[valU_yy; 0] - matK_2*[valU_yyy; valTHETA_y]...
            - matK_3*[valU_yyyy; valTHETA_yy]);
        
        vecDEFLECTION(valNSELE+3) = 2*vecDEFLECTION(valNSELE+2)...
            -vecDEFLECTION(valNSELE+1); % BC for deflection one element beyond wing (positive span direction)
        vecDEFLECTION(valNSELE+4) = 3*vecDEFLECTION(valNSELE+2)...
            -2*vecDEFLECTION(valNSELE+1); % BC for deflection two elements beyond wing (positive span direction)
        
        vecDEFLECTION(2) = vecDEFLECTION(4); % BC for deflection one element beyond root (negative span direction)
        vecDEFLECTION(1) = vecDEFLECTION(5); % BC for deflection two elements beyond root (negative span direction)
        
        vecTWIST(valNSELE+3) = vecTWIST(valNSELE+1); % BC for twist one element beyond wing (positive span direction)
        vecTWIST(2) = vecTWIST(4);
        vecTWIST(1) = vecTWIST(5);
        
        % Output result of deflection and twist to separate vectors
        vecDEFLECTION(yy) = tempTWISTBEND(1,:);
        vecTWIST(yy) = tempTWISTBEND(2,:);

    end

end

end