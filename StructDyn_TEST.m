
valSPAN = valSPAN/2;

valINITCOND = 0.0;

valDY = sum(2*vecDVEHVSPN,1)/length(vecDVEHVSPN);

valNSELE = sum(vecN,1)+1;

valSTRUCTDELTIME = valDELTIME;

C = -0.0333*vecSPANDIST + 0.76*ones(1,length(vecSPANDIST)) ;
tk = 0.02 ;
Tk = 0.15 ;
vecSPANAREA = pi*tk*C*(1 + Tk);

vecTWIST = zeros(1,valNSELE+4);

for valTIMESTEP = 3:3500
    
valSTRUCTTIME = valTIMESTEP;

if valTIMESTEP >= 3

    matDEF(1:valSTRUCTTIME-1,:) = matDEFGLOB(1:valTIMESTEP-1,:);

    matTWIST(1:valSTRUCTTIME-1,:) = matTWISTGLOB(1:valTIMESTEP-1,:);

    vecDEF(3) = 0; % Zero deflection at root BC
    vecTWIST(3) = 0; % Zero twist at root BC

    % Assemble load matrix
    matLOAD = [vecLIFTDIST' - vecLM.*9.81, vecMOMDIST - vecLM.*vecLSM.*9.81];

    for yy = 4:(valNSELE+2)

        %% Geometric property assembly

        % Assemble mass matrix
        matMASS = [vecLM(yy-2), -vecLM(yy-2).*vecLSM(yy-2); -vecLM(yy-2).*vecLSM(yy-2), vecLM(yy-2).*(vecLSM(yy-2)+(vecJT(yy-2)./(vecSPANAREA(yy-2))))];

        % Assemble stiffness matrices
        matK_1 = [matEIx(yy-2,3), 0; 0, 0];
        matK_2 = [matEIx(yy-2,2), 0; 0, -matGJt(yy-2,2)]; 
        matK_3 = [matEIx(yy-2,1), 0; 0, -matGJt(yy-2,1)];
        matB = [0 0; 0 10000];

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
        valTHETA_y = (matTWIST(valSTRUCTTIME-1,yy+1) - matTWIST(valSTRUCTTIME-1,yy-1))/2*valDY;
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

end

% Spanwise deflection and twist wrt structural timestep
vecDEF = matDEF(end,:);
vecTWIST = matTWIST(end,:);

% Spanwise deflection and twist wrt to global timestep

if valTIMESTEP >= 3
    
    matDEFGLOB(valTIMESTEP,:) = matDEF(end,:);
    matTWISTGLOB(valTIMESTEP,:) = matTWIST(end,:);

    matSLOPE(valTIMESTEP,:) = vecSLOPE';

else
    
    matDEFGLOB(valTIMESTEP,:) = matDEF(end,:);
    matTWISTGLOB(valTIMESTEP,:) = matTWIST(end,:);

    matSLOPE(valTIMESTEP,:) = vecSLOPE;

end

end

% end

