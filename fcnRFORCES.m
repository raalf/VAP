function [valCT, valCQ, valCP, vecCTCONV, vecCQCONV, vecCPCONV, vecDISTHRUST, vecDISNORM] = fcnRFORCES(valAZNUM, valDIA, valRPM, valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecQARM,vecDVEPITCH, vecDVETE, vecDVEWING, vecWDVEWING, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, vecCTCONV, vecCPCONV, vecCQCONV, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF, matTEPTS, matUINFTE)
% Force function where each force contributation is calculated and added.

%% Element normal forces, lift forces and side forces (freestream and induced)
[nind, nfree, thrustind, thrustfree, sideind, sidefree, axialind, axialfree,A,B,C] = fcnRDVENFORCE(valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEPITCH, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF);

% Calcualte thrust due to crossflow
[nfreecs,thrustCFfree, axialCFfree, sideCFfree] = fcnRCROSSFLOWFORCE(valNELE, vecTHETA, vecDVETESWP,vecDVELESWP,vecDVEHVSPN,vecDVEHVCRD,matDVE,matUINF,matVLST,B,C);

%% Calculate induced drag on the rotor
[inddrag, thrustinddrag, sideinddrag] = fcnRDVEINDDRAG(matCOEFF,matDVE,matVLST,matUINFTE,vecDVEHVSPN,vecDVETE, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN,vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP,vecSYM,vecDVEWING,vecWDVEWING, vecTHETA, matTEPTS);

%% Sum forces and calculate coefficients
[valCT, valCQ, valCP, vecCTCONV, vecCQCONV, vecCPCONV, vecDISTHRUST, vecDISNORM] = fcnROTORFORCE(nind, nfree, nfreecs, thrustind, thrustfree, thrustCFfree,  thrustinddrag, axialCFfree, axialind, axialfree, inddrag, sideind, sidefree, sideCFfree, sideinddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV,vecCQCONV, vecCPCONV, vecQARM, vecDVETE, vecTHETA);

% thrust = thrustind+thrustfree+thrustCFfree;
% torque = axialind+axialfree+inddrag;
% torque(18:34) = -1*torque(18:34);
% quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),0*thrust,0*thrust,thrust)
% quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),torque,0*torque,0*torque)
end 
