function [vecCLPDIST, vecCDPDIST, thrustind, thrustfree, thrustCFfree, tempTi, Pthrust, difthrustP, valCT, valFy, valFx, valCQ, valCP, valCMy, valCMx, vecCTCONV, vecCFyCONV, vecCFxCONV, vecCQCONV, vecCPCONV, vecCMyCONV, vecCMxCONV, vecDISTHRUST, vecDISNORM, vecDISAXIAL, vecDISSIDE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE, matWUINF, gamma_old, dGammadt] = fcnRFORCES(flagVERBOSE, valKINV, valMAXTIME, valAZNUM, valDIA, valRPM, valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEAREA, vecAIRFOIL, vecN, vecM, vecDVEPANEL, vecQARM,vecDVEPITCH, vecDVETE, vecDVEWING, vecWDVEWING, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, vecCTCONV, vecCFxCONV, vecCFyCONV, vecCPCONV, vecCMxCONV, vecCMyCONV, vecCQCONV, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF, matTEPTS, matUINFTE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE, flagSTEADY, gamma_old, dGammadt, flagVISCOUS, valNUMRO, vecDVEROTOR)
% Force function where each force contributation is calculated and added.

%% Element normal forces, lift forces and side forces (freestream and induced)
[nind, nfree, thrustind, thrustfree, sideind, sidefree, axialind, axialfree, A, B, C, matWUINF,  es, ea, en, gamma_old, dGammadt] = fcnRDVENFORCE(valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEPITCH, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF, flagSTEADY, gamma_old, dGammadt, valRPM, valAZNUM);

% Calcualte thrust due to crossflow
[nfreecs,thrustCFfree, axialCFfree, sideCFfree] = fcnRCROSSFLOWFORCE(valNELE, vecTHETA, vecDVETESWP,vecDVELESWP,vecDVEHVSPN,vecDVEHVCRD,matDVE,matUINF,matVLST,B,C);

%% Calculate induced drag on the rotor
[inddrag, thrustinddrag, sideinddrag] = fcnRDVEINDDRAG(matCOEFF,matDVE,matVLST,matUINFTE,vecDVEHVSPN,vecDVETE, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN,vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP,vecSYM,vecDVEWING,vecWDVEWING, vecTHETA, matTEPTS, flagSTEADY);

%% Apply viscous effects
if valTIMESTEP > valMAXTIME - valAZNUM && flagVISCOUS == 1% Only run for last rotation
    vecDISNORM = nfree + nind + nfreecs;
    induced_vel = fcnINDVEL(matCENTER,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
    vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
    matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
    vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, flagSTEADY);
    
    [vecCLPDIST, vecCDPDIST, Pthrust, Ptorque] = fcnRVISCOUS3(flagVERBOSE, ...
    valRPM,  valKINV, vecQARM, vecDVEHVCRD, vecN, vecM, ...
    vecDVELE, vecDVEPANEL, vecAIRFOIL, vecTHETA, vecDISNORM, vecDVEAREA,...
    matUINF, matVLST, matDVE, induced_vel);
%     [vecCLPDIST, vecCDPDIST, Pthrust, Ptorque, difthrustP, diffsideP, diffaxialP] = fcnRVISCOUS(en, es, ea,  flagVERBOSE, ...
%         valRPM, valDIA, valKINV, vecQARM, vecDVEHVCRD, vecN, vecM, ...
%         vecDVELE, vecDVEPANEL, vecAIRFOIL, vecTHETA, vecDISNORM, vecDVEAREA,...
%         matUINF, matVLST, matDVE, induced_vel);

%     [difthrustP, diffaxialP] = fcnRVISCOUS2(flagVERBOSE, valRPM, valKINV, vecQARM, vecTHETA, vecDVEAREA, vecAIRFOIL, vecDVEHVCRD, vecN, vecM, vecDISNORM, vecDVELE, vecDVEPANEL, matVLST, matDVE, matUINF, induced_vel);
%     vecCDPDIST = [];
%     vecCLPDIST = [];
%     Pthrust = zeros(valNELE,1);
%     Ptorque = zeros(valNELE,1);
    difthrustP = zeros(valNELE,1);
    diffsideP = zeros(valNELE,1);
    diffaxialP = zeros(valNELE,1);
else
    Pthrust = zeros(valNELE,1);
    Ptorque = zeros(valNELE,1);
    difthrustP = zeros(valNELE,1);
    diffsideP = zeros(valNELE,1);
    diffaxialP = zeros(valNELE,1);
    vecCDPDIST = [];
    vecCLPDIST = [];
end
%% Sum forces and calculate coefficients
[thrustind, thrustfree, thrustCFfree, tempTi, Pthrust, difthrustP, valCT, valFy, valFx, valCQ, valCP, valCMy, valCMx, vecCTCONV, vecCFyCONV, vecCFxCONV, vecCQCONV, vecCPCONV, vecCMyCONV, vecCMxCONV, vecDISTHRUST, vecDISNORM, vecDISAXIAL, vecDISSIDE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE] = fcnROTORFORCE(difthrustP, diffsideP, diffaxialP, nind, nfree, nfreecs, Pthrust, thrustind, thrustfree, thrustCFfree,  thrustinddrag, Ptorque, axialCFfree, axialind, axialfree, inddrag, sideind, sidefree, sideCFfree, sideinddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV, vecCFxCONV, vecCFyCONV, vecCQCONV, vecCPCONV, vecCMxCONV, vecCMyCONV, vecQARM, vecDVETE, vecTHETA, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE, valNUMRO, vecDVEROTOR);
% thrust = thrustind+thrustfree+thrustCFfree;
% torque = axialind+axialfree+inddrag;
% torque(18:34) = -1*torque(18:34);
% quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),0*thrust,0*thrust,thrust)
% quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),torque,0*torque,0*torque)
end 
