function [THRUSTDIST, TORQUEDIST] = fcnRVISCOUS2(flagVERBOSE, valRPM, valKINV, vecQARM, vecTHETA, vecDVEAREA, vecAIRFOIL, vecDVEHVCRD, vecN, vecM, vecDISNORM, vecDVELE, vecDVEPANEL, matVLST, matDVE, matUINF, matINDVEL)
% This is an updated viscous correction for rotors. This function uses
% effective angle of attack to search though the lookup tables instead of
% the Cl values

% Calculate chordline direction at midspan of each dve
avgle = (matVLST(matDVE(:,1),:)+matVLST(matDVE(:,2),:))./2;
avgte = (matVLST(matDVE(:,3),:)+matVLST(matDVE(:,4),:))./2;
tempDIF = avgte - avgle;
matCRDLINE = (tempDIF)./(repmat(sqrt(sum(tempDIF.^2,2)),[1,3]));
% Calculate velocity 
vecV = dot(matUINF + matINDVEL, matCRDLINE,2);

% Calculate effective angle of attack
tempUINF = matUINF +  matINDVEL;
tempVELMAG = sqrt(tempUINF(:,1).^2 + tempUINF(:,2).^2 + tempUINF(:,3).^2);
tempCRDMAG = sqrt(matCRDLINE(:,1).^2 + matCRDLINE(:,2).^2 + matCRDLINE(:,3).^2);
vecALPHAEFF = acos(dot(tempUINF, matCRDLINE,2)./(tempVELMAG.*tempCRDMAG));

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

vecREDIST = vecV.*2.*sum(vecDVEHVCRD(rows),2)./valKINV;
len = 0;

vecCDPDIST = zeros(size(rows,1),1);
vecCLPDIST = zeros(size(rows,1),1);
for j = 1:length(idxpanel)
    pan = idxpanel(j);
    airfoil = dlmread(strcat('airfoils/airfoil',num2str(vecAIRFOIL(pan)),'.dat'),'', 1, 0);

    HiRe = airfoil(end,4);
    LoRe = airfoil(1,4);
    
    alpha = vecALPHAEFF(len + j);
    %cl = vecCNDIST(len + j);

    if vecREDIST(len + j) > HiRe
        if flagVERBOSE == 1
            fprintf('\nRe higher than airfoil Re data')
        end
        Re2 = airfoil(end,4);
        temp_var = airfoil(airfoil(:,4) == Re2, 1);
        cl_max_alpha = temp_var(end);
    elseif vecREDIST(len + j) < LoRe
        if flagVERBOSE == 1
            fprintf('\nRe lower than airfoil Re data');
        end
        Re2 = airfoil(1,4);
        temp_var = airfoil(airfoil(:,4) == Re2, 1);
        cl_max_alpha = temp_var(end);
    else
        re1 = airfoil(airfoil(:,4) < vecREDIST(len + j), 4);
        re1 = re1(end);
        cl_max_alpha1 = airfoil(airfoil(:,4) < vecREDIST(len + j), 1);
        cl_max_alpha1 = cl_max_alpha1(end);

        temp_var = airfoil(airfoil(:,4) > vecREDIST(len + j),4);
        re2 = temp_var(1);
        temp_var = airfoil(airfoil(:,4) == (temp_var(1)), 1);
        cl_max_alpha2 = temp_var(end);

        cl_max_alpha = interp1([re1 re2],[cl_max_alpha1 cl_max_alpha2], vecREDIST(len + j));
    end
    % correcting the section cl if we are above cl_max
    if radtodeg(alpha) > cl_max_alpha
        if flagVERBOSE == 1
            fprintf('\nBlade Stall on Section %d, alpha = %f Re = %0.0f', j, radtodeg(alpha), vecREDIST(len + j))
        end
        %vecCNDIST0(len+j) = 0.825*cl_max; % setting the stalled 2d cl     
        
        % Make apply stall model using empirical equations
        % cn = cd,90*(sin(alpha_eff))/(0.56+0.44sin(alpha_eff))
        % ct = cd,0*cos(alpha_eff)/2
        % cd = cn*sin(alpha_eff)+ct*cos(alpha_eff)
        % Note: cd,90 = 2
        % Find cd_0
        temp = scatteredInterpolant(airfoil(:,4), airfoil(:,1), airfoil(:,3),'nearest');
        cd_0 = temp(vecREDIST(len+j),0);
        %alpha_eff = vecCNDIST/(2*pi);
        %alpha_eff = asin((vecCNDIST(len+j)/2)*0.56/(1-0.44*((vecCNDIST(len+j)/2))));
        cn = 2*sin(abs(vecALPHAEFF(len+j)))/(0.56+0.44*sin(abs(vecALPHAEFF(len+j))));
        ct = cd_0*cos(abs(vecALPHAEFF(len+j)))/2;
        
        vecCDPDIST(len + j) = cn*sin(abs(vecALPHAEFF(len+j))) + ct*cos(abs(vecALPHAEFF(len+j)));
        vecCLPDIST(len + j) = cn*cos(abs(vecALPHAEFF(len+j))) - ct*sin(abs(vecALPHAEFF(len+j)));
    else
        warning off
        F = scatteredInterpolant(airfoil(:,4), airfoil(:,1), airfoil(:,3),'nearest');
        vecCDPDIST(len + j, 1) = F(vecREDIST(len + j), radtodeg(alpha));
    end
end
% Calculate viscous drag distribution
vecDPDIST = 0.5*(vecCDPDIST.*((vecV).^2).*(sum(vecDVEAREA(rows),2)));

% Apply an appropriate directions
tempUINFX = tempUINF(:,1);
tempUINFY = tempUINF(:,2);
tempUINFZ = tempUINF(:,3);
matUINFAVG = [sum(tempUINFX(rows),2)/(size(rows,2)) sum(tempUINFY(rows),2)/(size(rows,2)) sum(tempUINFZ(rows),2)/(size(rows,2))];
tempDIR = [matUINFAVG(:,1).*cos(vecTHETA(rows(:,1))) matUINFAVG(:,2).*sin(vecTHETA(rows(:,1))) matUINFAVG(:,3)];
tempDIR = tempDIR./(sqrt(sum(tempDIR.^2,2)));
matDPDIST = vecDPDIST.*tempDIR;

% Resolve to viscous thrust and torque
THRUSTDIST = matDPDIST(:,3);
TORQUEDIST = (dot(matDPDIST,[abs(cos(vecTHETA(rows(:,1)))) abs(sin(vecTHETA(rows(:,1)))) zeros(size(matDPDIST,1),1)],2)).*vecQARM(vecDVELE==1);
POWERDIST = TORQUEDIST*2.*pi.*(valRPM./60);

end