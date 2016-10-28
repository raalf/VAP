clear,clc

valTIMESTEP = 0;
valNSELE = 200;
vecLIFT = 5000*ones(1,valNSELE); 
vecMOM = zeros(1,valNSELE);
matDEFLECTION = zeros(1,valNSELE+4);
matTWIST = zeros(1,valNSELE+4);
vecSPANAREA = 1750*(1/1000)^2*ones(1,valNSELE);
vecSPANINERTIA = (1/1000)^4*7414583.333333332*ones(1,valNSELE);
vecSPANINERTIA_1 = zeros(1,valNSELE);
vecSPANINERTIA_2 = zeros(1,valNSELE);
matSPANINERTIA = [vecSPANINERTIA; vecSPANINERTIA_1; vecSPANINERTIA_2];
vecLINMASS = 2700*vecSPANAREA; 
vecTORSIONRIGIDITY = 27e9*14583.3333333333*(1/1000)^4*ones(1,valNSELE);
vecMASS2SHEAR = 0.1*ones(1,valNSELE);
valYMODULUS = 69e9; 
valSPAN = 5;

valTOL = 100;

valINITCOND = 0.0;

% Preallocate space for deflection and twist vectors
vecDEFLECTION = zeros(1,valNSELE+4);
vecTWIST = zeros(1,valNSELE+4);

% Span increment between each element (in m)
valDY = valSPAN/(valNSELE-1);

% Vector of beam spanwise locations where deflection and twist is computed
vecBEAM = 0:valDY:valSPAN;

valDELTIME = sqrt(150000*valDY^4/valYMODULUS);

% while valTOL > 0.00005

for valTIMESTEP = 1:5000
    
    % valTIMESTEP = valTIMESTEP + 1;

    [vecDEFLECTION, vecTWIST, vecBEAM] = fcnWINGTWISTBEND(valDELTIME, valTIMESTEP, vecLIFT, vecMOM, matDEFLECTION, matTWIST, vecSPANAREA, matSPANINERTIA, vecLINMASS, vecTORSIONRIGIDITY, vecMASS2SHEAR, valYMODULUS, valNSELE, valSPAN, vecDEFLECTION, vecTWIST, vecBEAM, valINITCOND, valDY);
    
    matDEFLECTION(valTIMESTEP,:) = vecDEFLECTION;
    matTWIST(valTIMESTEP,:) = vecTWIST;

    % if valTIMESTEP > 2
	   %  valTOL = abs((matDEFLECTION(valTIMESTEP,valNSELE+2)-matDEFLECTION(valTIMESTEP-1,valNSELE+2))/...
	   %  	matDEFLECTION(valTIMESTEP-1,valNSELE+2))*100;
    % end
    
end

figure(1)
plot(vecBEAM(:),matDEFLECTION(valTIMESTEP,3:valNSELE+2))
xlabel('Span (m)')
ylabel('Deflection (m)')
title('Wing Spanwise Deflection')

figure(2)
plot(vecBEAM(:),(180/(2*pi)).*matTWIST(valTIMESTEP,3:valNSELE+2))
xlabel('Span (m)')
ylabel('Twist (deg)')
title('Wing Spanwise Twist')