function [valCT, valCQ, valCP] = fcnRVISCOUS(flagVERBOSE, valCT, valCQ, ...
    valCP, valRPM, valDIA, valKINV, vecQARM, vecDVEHVCRD, vecN, vecM, ...
    vecDVELE, vecDVEPANEL, vecAIRFOIL, vecTHETA, vecDISNORM, vecDVEAREA,...
    matUINF, matVLST, matDVE)
% This function applies a viscous correction using look up tables.
% OUTPUT
%   valCT - Viscous corrected thrust coeff
%   valCQ - Viscous corrected torque coeff
%   valCP - Viscous corrected power coeff


% Note: CN is normal to both the freestream and the spanwise direction
% Calculate velocity seen by section
% OLD velocity calculation
% tempXVEL = matUINF(:,1).*cos(vecTHETA)+matUINF(:,2).*sin(vecTHETA);
% vecV1 = sqrt(tempXVEL.^2 + matUINF(:,3).^2);

% Calculate chordline direction at midspan of each dve
avgle = (matVLST(matDVE(:,1),:)+matVLST(matDVE(:,2),:))./2;
avgte = (matVLST(matDVE(:,3),:)+matVLST(matDVE(:,4),:))./2;
tempDIF = avgte - avgle;
matCRDLINE = (tempDIF)./(repmat(sqrt(sum(tempDIF.^2,2)),[1,3]));
% Calculate velocity 
vecV = dot(matUINF, matCRDLINE,2);

[ledves, ~, ~] = find(vecDVELE > 0);
lepanels = vecDVEPANEL(ledves);

% vecDVEROTOR is which rotor each DVE is associated to and must be updated
% for mutiple rotors
len = size(vecDISNORM,1);
vecDVEROTOR = ones(len,1);
idxdve = ledves(vecDVEROTOR(ledves) == 1);
idxpanel = lepanels(vecDVEROTOR(ledves) == 1);

m = vecM(idxpanel);

% Matrix of how much to add to an index to get the next chordwise element
tempm = repmat(vecN(idxpanel),1, m(1)).*repmat([0:m(1)-1],length(idxpanel),1);

% Which row index to add
rows = repmat(idxdve,1,m(1)) + tempm;

% Average velocities across chord
vecV = sum(vecV(rows),2)/(size(rows,2));

% CN = 2*(N/rho)/(V^2*S)
vecCNDIST = (sum(vecDISNORM(rows),2).*2)./(sum(vecDVEAREA(rows),2).*(vecV.^2));
vecREDIST = vecV.*2.*sum(vecDVEHVCRD(rows),2)./valKINV;

vecCDPDIST = zeros(size(rows,1),1);
len = 0;
for j = 1:length(idxpanel)
    pan = idxpanel(j);
    airfoil = dlmread(strcat('airfoils/airfoil',num2str(vecAIRFOIL(pan)),'.dat'),'', 1, 0);

    HiRe = airfoil(end,4);
    LoRe = airfoil(1,4);

    cl = vecCNDIST(len + j);

    if vecREDIST(len + j) > HiRe
        if flagVERBOSE == 1
            fprintf('\nRe higher than airfoil Re data')
        end
        Re2 = airfoil(end,4);
        temp_var = airfoil(airfoil(:,4) == Re2, 2);
        cl_max = temp_var(end);
    elseif vecREDIST(len + j) < LoRe
        if flagVERBOSE == 1
            fprintf('\nRe lower than airfoil Re data');
        end
        Re2 = airfoil(1,4);
        temp_var = airfoil(airfoil(:,4) == Re2, 2);
        cl_max = temp_var(end);
    else
        re1 = airfoil(airfoil(:,4) < vecREDIST(len + j), 4);
        re1 = re1(end);
        cl_max1 = airfoil(airfoil(:,4) < vecREDIST(len + j), 2);
        cl_max1 = cl_max1(end);

        temp_var = airfoil(airfoil(:,4) > vecREDIST(len + j),4);
        re2 = temp_var(1);
        temp_var = airfoil(airfoil(:,4) == (temp_var(1)), 2);
        cl_max2 = temp_var(end);

        cl_max = interp1([re1 re2],[cl_max1 cl_max2], vecREDIST(len + j));
    end

    % correcting the section cl if we are above cl_max
    if cl > cl_max
        if flagVERBOSE == 1
            fprintf('\nBlade Stall on Section %d, cl = %f Re = %0.0f', 1, j, cl, vecREDIST(len + j))
        end
        vecCNDIST(len + j) = 0.825*cl_max; % setting the stalled 2d cl
    end

    F = scatteredInterpolant(airfoil(:,4), airfoil(:,2), airfoil(:,3),'nearest');
    vecCDPDIST(len + j, 1) = F(vecREDIST(len + j), cl);
end
% Calculate viscous drag distribution
vecDPDIST = 0.5*(vecCDPDIST.*((vecV).^2).*(sum(vecDVEAREA(rows),2)));

% Apply an appropriate directions
tempUINFX = matUINF(:,1);
tempUINFY = matUINF(:,2);
tempUINFZ = matUINF(:,3);
matUINFAVG = [sum(tempUINFX(rows),2)/(size(rows,2)) sum(tempUINFY(rows),2)/(size(rows,2)) sum(tempUINFZ(rows),2)/(size(rows,2))];
tempDIR = [matUINFAVG(:,1).*cos(vecTHETA(rows(:,1))) matUINFAVG(:,2).*sin(vecTHETA(rows(:,1))) matUINFAVG(:,3)];
tempDIR = tempDIR./(sqrt(sum(tempDIR.^2,2)));
matDPDIST = vecDPDIST.*tempDIR;

% Resolve to viscous thrust and torque
THRUSTDIST = matDPDIST(:,3);
TORQUEDIST = (dot(matDPDIST,[abs(cos(vecTHETA(rows(:,1)))) abs(sin(vecTHETA(rows(:,1)))) zeros(size(matDPDIST,1),1)],2)).*vecQARM(vecDVELE==1);
POWERDIST = TORQUEDIST*2.*pi.*(valRPM./60);

% Calculate coefficient
valCTP = (sum(THRUSTDIST))/(((valRPM/60)^2)*((valDIA)^4));
valCQP = (sum(TORQUEDIST))/(((valRPM/60)^2)*((valDIA)^5));
valCPP = (sum(POWERDIST))/((valRPM/60)^3*(valDIA^5));

% Add viscous forces to thrust and power
valCT = valCT + valCTP;
valCQ = valCQP + valCQ;
valCP = valCPP + valCP;
end