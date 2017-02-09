function [nind, nfree] = fcnRFORCES( valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEPITCH, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF)
% Force function where each force contributation is calculated and added.

%% Element normal forces, lift forces and side forces (freestream and induced)
[nind, nfree, thrustind, thrustfree, sideind, sidefree, axialind, axialfree] = fcnRDVENFORCE(valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEPITCH, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF);

%% Calculate induced drag on the rotor


end

