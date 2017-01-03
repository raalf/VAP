function [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF] = fcnWINGTWISTBEND(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, matGJt, vecLSM, vecN, valSPAN, vecDVEHVSPN, valTIMESTEP, matDEFGLOB, matTWISTGLOB)
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

valSPAN = valSPAN/2;

valINITCOND = 0.0;

valDY = (valSPAN/2)/sum(vecN,1);

valNSELE = sum(vecN,1)+1;

valDELTIME = sqrt(1.5*valDY^4/710^2);

vecBEAM = 0:valDY:valSPAN/2;

Jy = (9e-5*vecBEAM - 0.0043*vecBEAM + 0.0486)' ;

vecDEF = zeros(1,valNSELE+4);
vecTWIST = zeros(1,valNSELE+4);

vecLIFTDIST = [vecLIFTDIST; 0];

vecMOMDIST = [vecMOMDIST; 0];

% Temporary area calculation
C = -0.0333*vecBEAM + 0.76*ones(1,length(vecBEAM)) ;
tk = 0.02 ;
Tk = 0.15 ;
vecSPANAREA = pi*tk*C*(1 + Tk) ;

% Temporary AC location calculation
LSAC = 0.0062*vecBEAM.*vecBEAM.*vecBEAM - 0.0533*vecBEAM.*vecBEAM + 0.1403*vecBEAM + 0.7029;

vecLSM = 0.1*LSAC;

%% Beam boundary conditions

for valSTRUCTTIME = 1:30000
   
    if valSTRUCTTIME == 1
        
        % Initial conditions for wing deflection
        vecDEF(3:valNSELE+2) = valINITCOND*vecBEAM.*vecBEAM;
        vecDEF(2) = vecDEF(4); % BC for deflection one element beyond root (negative span direction)

        vecDEF(valNSELE+3) = 2*vecDEF(valNSELE+2)...
            -vecDEF(valNSELE+1); % BC for deflection one element beyond wing (positive span direction)
        vecDEF(valNSELE+4) = 3*vecDEF(valNSELE+2)...
            -2*vecDEF(valNSELE+1); % BC for deflection two elements beyond wing (positive span direction)   
        
        % Initial conditions for wing twist
        vecTWIST(3:valNSELE+2) = valINITCOND*(vecBEAM./valSPAN);
        vecTWIST(valNSELE+3) = vecTWIST(valNSELE+1); % BC for twist one element beyond wing (positive span direction)

        matDEF(valSTRUCTTIME,:) = vecDEF;
        matTWIST(valSTRUCTTIME,:) = vecTWIST;

    elseif valSTRUCTTIME == 2
        
        % Initial conditions for wing deflection
        vecDEF(3:valNSELE+2) = valINITCOND*vecBEAM.*vecBEAM;
        vecDEF(2) = vecDEF(4); % BC for deflection one element beyond root (negative span direction)
        
        vecDEF(valNSELE+3) = 2*vecDEF(valNSELE+2)...
            -vecDEF(valNSELE+1); % BC for deflection one element beyond wing (positive span direction)
        vecDEF(valNSELE+4) = 3*vecDEF(valNSELE+2)...
            -2*vecDEF(valNSELE+1); % BC for deflection two elements beyond wing (positive span direction)
        
        % Initial conditions for wing twist
        vecTWIST(3:valNSELE+2) = valINITCOND*(vecBEAM./valSPAN);
        vecTWIST(valNSELE+3) = vecTWIST(valNSELE+1); % BC for twist one element beyond wing (positive span direction)

        matDEF(valSTRUCTTIME,:) = vecDEF;
        matTWIST(valSTRUCTTIME,:) = vecTWIST;

    else
                
        vecDEF(3) = 0; % Zero deflection at root BC
        vecTWIST(3) = 0; % Zero twist at root BC
        
        for yy = 4:(valNSELE+2)
            
            %% Geometric property assembly

            % Assemble mass matrix
            matMASS = [vecLM(yy-2), -vecLM(yy-2).*vecLSM(yy-2); -vecLM(yy-2).*vecLSM(yy-2), vecLM(yy-2).*(vecLSM(yy-2)+(Jy(yy-2)./(vecSPANAREA(yy-2))))];

            % Assemble load matrix
            matLOAD = [vecLIFTDIST(yy-2); vecMOMDIST(yy-2)];

            % Assemble stiffness matrices
            matK_1 = [matEIx(yy-2,3), 0; 0, 0];
            matK_2 = [matEIx(yy-2,2), 0; 0, -matGJt(yy-2,2)]; 
            matK_3 = [matEIx(yy-2,1), 0; 0, -matGJt(yy-2,1)];
            matB = [0 0;0 100];

            %% Finite difference relations for partial derivatives

            % Finite difference relations for partial derivative of deflection w.r.t Y
            valU_yy = (matDEF(valSTRUCTTIME-1,yy+1) - 2*matDEF(valSTRUCTTIME-1,yy) + matDEF(valSTRUCTTIME-1,yy-1))/(valDY)^2;
            valU_yyy = (matDEF(valSTRUCTTIME-1,yy+2) - 3*matDEF(valSTRUCTTIME-1,yy+1) + 3*matDEF(valSTRUCTTIME-1,yy)- ...
                matDEF(valSTRUCTTIME-1,yy-1))/(valDY)^3;
            valU_yyyy = (matDEF(valSTRUCTTIME-1,yy+2) - 4*matDEF(valSTRUCTTIME-1,yy+1) + 6*matDEF(valSTRUCTTIME-1,yy) - ...
                4*matDEF(valSTRUCTTIME-1,yy-1) + matDEF(valSTRUCTTIME-1,yy-2))/(valDY)^4;
            
            valTDOT = (matTWIST(valSTRUCTTIME-1,yy) - matTWIST(valSTRUCTTIME - 2,yy))./valDELTIME;
            valHDOT = (matDEF(valSTRUCTTIME-1,yy) - matDEF(valSTRUCTTIME - 2, yy))./valDELTIME;

            % Finite difference relations for partial derivative of twist w.r.t Y
            valTHETA_y = (matTWIST(valSTRUCTTIME-1,yy+1) - matTWIST(valSTRUCTTIME-1,yy-1))/2*valDY;
            valTHETA_yy = (matTWIST(valSTRUCTTIME-1,yy+1) - 2*matTWIST(valSTRUCTTIME-1,yy) + matTWIST(valSTRUCTTIME-1,yy-1))/(valDY^2);

            %% Solve matrix equation

            % Temp variable with the wing deflection and twist stored as a matrix. The
            % first row is the deflection, w/ each column as a spanwise station. The
            % second row is the twist, w/ each column as a spanwise station.

            tempTWISTBEND = 2.*[matDEF(valSTRUCTTIME-1,yy); matTWIST(valSTRUCTTIME-1,yy)] - [matDEF(valSTRUCTTIME-2,yy); matTWIST(valSTRUCTTIME-2,yy)] ...
                + (valDELTIME^2).*inv(matMASS)*(matLOAD - matK_1*[valU_yy; 0] - matK_2*[valU_yyy; valTHETA_y] -...
                matK_3*[valU_yyyy; valTHETA_yy]);

            % Output result of deflection and twist to separate vectors
            vecDEF(yy) = tempTWISTBEND(1,:);
            vecTWIST(yy) = tempTWISTBEND(2,:);

        end
        
        vecDEF(valNSELE+3) = 2*vecDEF(valNSELE+2)...
            -vecDEF(valNSELE+1); % BC for deflection one element beyond wing (positive span direction)
        vecDEF(valNSELE+4) = 3*vecDEF(valNSELE+2)...
            -2*vecDEF(valNSELE+1); % BC for deflection two elements beyond wing (positive span direction)

        vecDEF(2) = vecDEF(4); % BC for deflection one element beyond root (negative span direction)
            
        vecTWIST(valNSELE+3) = vecTWIST(valNSELE+1); % BC for twist one element beyond wing tip (positive span direction)

    end
    
    matDEF(valSTRUCTTIME,:) = vecDEF;
    matTWIST(valSTRUCTTIME,:) = vecTWIST;

end

% Spanwise deflection and twist wrt structural timestep
vecDEF = matDEF(valSTRUCTTIME,:);
vecTWIST = matTWIST(valSTRUCTTIME,:);

% Spanwise deflection and twist wrt to global timestep
matDEFGLOB(valTIMESTEP,:) = vecDEF;
matTWISTGLOB(valTIMESTEP,:) = vecTWIST;

end