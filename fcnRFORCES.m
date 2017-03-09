function [valCT, valCQ,vecCTCONV] = fcnRFORCES(valAZNUM, valDIA, valRPM, valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecCPRADI,vecDVEPITCH, vecDVETE, vecDVEWING, vecWDVEWING, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, vecCTCONV, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF, matTEPTS, matUINFTE)
% Force function where each force contributation is calculated and added.

%% Element normal forces, lift forces and side forces (freestream and induced)
[nind, nfree, thrustind, thrustfree, sideind, sidefree, axialind, axialfree,A,B,C] = fcnRDVENFORCE(valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEPITCH, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF);

% Calcualte thrust due to crossflow
[nfreecs,thrustCFfree] = fcnRCROSSFLOWFORCE(vecDVETESWP,vecDVELESWP,vecDVEHVSPN,vecDVEHVCRD,matDVE,matUINF,matVLST,B,C);

%% Calculate induced drag on the rotor
[inddrag] = fcnRDVEINDDRAG(matCOEFF,matDVE,matVLST,matUINFTE,vecDVEHVSPN,vecDVETE, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN,vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP,vecSYM,vecDVEWING,vecWDVEWING, vecTHETA, matTEPTS);

%% Sum forces and calculate coefficients
[valCT, valCQ, vecCTCONV] = fcnROTORFORCE(thrustind, thrustfree, thrustCFfree, axialind, axialfree, inddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV, vecCPRADI);
end 
