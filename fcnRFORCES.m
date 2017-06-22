function [valCT, valFy, valFx, valCQ, valCP, valCMy, valCMx, vecCTCONV, vecCFyCONV, vecCFxCONV, vecCQCONV, vecCPCONV, vecCMyCONV, vecCMxCONV, vecDISTHRUST, vecDISNORM, vecDISAXIAL, vecDISSIDE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE, matWUINF] = fcnRFORCES(flagVERBOSE, valKINV, valMAXTIME, valAZNUM, valDIA, valRPM, valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEAREA, vecAIRFOIL, vecN, vecM, vecDVEPANEL, vecQARM,vecDVEPITCH, vecDVETE, vecDVEWING, vecWDVEWING, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, vecCTCONV, vecCFxCONV, vecCFyCONV, vecCPCONV, vecCMxCONV, vecCMyCONV, vecCQCONV, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF, matTEPTS, matUINFTE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE)
% Force function where each force contributation is calculated and added.

%% Element normal forces, lift forces and side forces (freestream and induced)
[nind, nfree, thrustind, thrustfree, sideind, sidefree, axialind, axialfree,A,B,C, matWUINF,  es, ea, en] = fcnRDVENFORCE(valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEPITCH, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF);

% Calcualte thrust due to crossflow
[nfreecs,thrustCFfree, axialCFfree, sideCFfree] = fcnRCROSSFLOWFORCE(valNELE, vecTHETA, vecDVETESWP,vecDVELESWP,vecDVEHVSPN,vecDVEHVCRD,matDVE,matUINF,matVLST,B,C);

%% Calculate induced drag on the rotor
[inddrag, thrustinddrag, sideinddrag] = fcnRDVEINDDRAG(matCOEFF,matDVE,matVLST,matUINFTE,vecDVEHVSPN,vecDVETE, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN,vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP,vecSYM,vecDVEWING,vecWDVEWING, vecTHETA, matTEPTS);

%% Apply viscous effects
if valTIMESTEP > valMAXTIME - valAZNUM % Only run for last rotation
    vecDISNORM = nfree + nind + nfreecs;
    
    [Pthrust, Ptorque, difthrustP, diffsideP, diffaxialP] = fcnRVISCOUS(en, es, ea,  flagVERBOSE, ...
        valRPM, valDIA, valKINV, vecQARM, vecDVEHVCRD, vecN, vecM, ...
        vecDVELE, vecDVEPANEL, vecAIRFOIL, vecTHETA, vecDISNORM, vecDVEAREA,...
        matUINF, matVLST, matDVE, matWUINF);
else
    Pthrust = zeros(valNELE,1);
    Ptorque = zeros(valNELE,1);
    difthrustP = zeros(valNELE,1);
    diffsideP = zeros(valNELE,1);
    diffaxialP = zeros(valNELE,1);
end
%% Sum forces and calculate coefficients
[valCT, valFy, valFx, valCQ, valCP, valCMy, valCMx, vecCTCONV, vecCFyCONV, vecCFxCONV, vecCQCONV, vecCPCONV, vecCMyCONV, vecCMxCONV, vecDISTHRUST, vecDISNORM, vecDISAXIAL, vecDISSIDE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE] = fcnROTORFORCE(difthrustP, diffsideP, diffaxialP, nind, nfree, nfreecs, Pthrust, thrustind, thrustfree, thrustCFfree,  thrustinddrag, Ptorque, axialCFfree, axialind, axialfree, inddrag, sideind, sidefree, sideCFfree, sideinddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV, vecCFxCONV, vecCFyCONV, vecCQCONV, vecCPCONV, vecCMxCONV, vecCMyCONV, vecQARM, vecDVETE, vecTHETA, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE);


% thrust = thrustind+thrustfree+thrustCFfree;
% torque = axialind+axialfree+inddrag;
% torque(18:34) = -1*torque(18:34);
% quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),0*thrust,0*thrust,thrust)
% quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),torque,0*torque,0*torque)
end 
