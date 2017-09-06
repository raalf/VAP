function [valDELTIME, matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST,...
    vecSPANDIST, matSC, vecMAC, vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF_old, matTWIST_old, matSLOPE,...
    matNPVLST, matNPNEWWAKE, matNEWWAKE, valUINF, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, matVLST, matDVE, matCENTER, matUINF, valGUSTTIME, flagSTEADY] = fcnFLEXWING(vecDVEHVSPN,...
    vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF, vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL,...
    vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecLIFTDIST, vecMOMDIST, valSPAN, valTIMESTEP, matDEFGLOB, matTWISTGLOB,...
    matSLOPE, valALPHA, valBETA, matVLST, matCENTER, matDVE, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF,...
    flagSTATIC, valSDELTIME, valDELTIME, matDEF, matTWIST, valSTIFFSTEPS, valGUSTTIME, valGUSTAMP, valGUSTL,...
    flagGUSTMODE, flagSTEADY, vecDVEHVCRD, vecLEDVES, vecLAMBDA)

% valDELTIME = valSDELTIME;

[matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST, matSC, vecMAC] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
    vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW, matCENTER, valSPAN);

matCENTER_old = matCENTER;

% Applies gust after aeroelastic convergence
% Michael A. D. Melville, Denver, CO, 80218
% if valGUSTTIME > 1
    
    valDELTIME = valSDELTIME;
    flagSTEADY = 2;
    valSTIMESTEP = 1;
    valDELTIME = 0.2;
    valMAXTIME = 1000000;
    
    matDEF = zeros(valMAXTIME,size(vecLM,1)-1);
    matTWIST = zeros(valMAXTIME,size(vecLM,1)-1);
    
    for valTIMESTEP = 4:valMAXTIME
    [matDEF, matTWIST] = fcnIMPLICIT(matEIx, matGJt, matDEF, matTWIST, vecJT, valDELTIME, vecDVEHVSPN, vecLSM, vecLM, vecLIFTDIST,...
        vecMOMDIST, vecSPANDIST, valTIMESTEP);
%     if valDELTIME > 0.01
%         valDELTIME = valDELTIME - 0.00001;
%     else
%         valDELTIME = 0.01;
%     end
    
    end       
    
    matDEF_old = matDEF;
    matTWIST_old = matTWIST;
    
[matNPVLST, matNPNEWWAKE, matNEWWAKE, valUINF] = fcnMOVEFLEXWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, vecDVEHVSPN, vecDVELE, matNPVLST, matDEFGLOB,...
    matTWISTGLOB, matSLOPE, valTIMESTEP, vecN, vecM, vecDVEWING, vecDVEPANEL, matSCLST, vecDVEPITCH, matNPDVE, vecSPANDIST, vecCL, valWEIGHT, valAREA, valDENSITY,valUINF);

[ vecDVEHVSPN, vecDVEHVCRD, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEAREA, matDVENORM, matVLST, matDVE, matCENTER, matNEWWAKE ] = fcnVLST2DVEPARAM( matNPDVE, matNPVLST, matNEWWAKE, vecDVETE );

[matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME, valTIMESTEP);

% Determine % relative change between aeroelastic timesteps
tol_def = (100*abs(matDEFGLOB(valTIMESTEP,end-2)-matDEFGLOB(valTIMESTEP-valSTIFFSTEPS,end-2))/matDEFGLOB(valTIMESTEP-valSTIFFSTEPS,end-2));
tol_twist = (100*abs(matTWISTGLOB(valTIMESTEP,end-2)-matTWISTGLOB(valTIMESTEP-valSTIFFSTEPS,end-2))/matTWISTGLOB(valTIMESTEP-valSTIFFSTEPS,end-2));

% Add in gust velocities to matUINF if convergence tolerance is met
if (tol_def < 0.05 && tol_twist < 0.05) || valGUSTTIME > 1
    [matUINF] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF);
    valGUSTTIME = valGUSTTIME + 1;
end

end