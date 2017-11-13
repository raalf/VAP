function [valDELTIME, matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST,...
    vecSPANDIST, matSC, vecMAC, vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF_old, matTWIST_old, matSLOPE,...
    matNPVLST, matNPNEWWAKE, matNEWWAKE, valUINF, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, matVLST, matDVE, matCENTER, matUINF, valGUSTTIME, flagSTEADY, gust_vel_old,zvel] = fcnFLEXWING(vecDVEHVSPN,...
    vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF, vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL,...
    vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecLIFTDIST, vecMOMDIST, valSPAN, valTIMESTEP, matDEFGLOB, matTWISTGLOB,...
    matSLOPE, valALPHA, valBETA, matVLST, matCENTER, matDVE, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF,...
    flagSTATIC, valSDELTIME, valDELTIME, valNSELE, matDEF, matTWIST, valSTIFFSTEPS, valGUSTTIME, valGUSTAMP, valGUSTL,...
    valGUSTSTART, flagGUSTMODE, flagSTEADY, gust_vel_old, matUINF, zvel)

% valDELTIME = valSDELTIME;

[matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST, matSC, vecMAC] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
    vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW);

matCENTER_old = matCENTER;

% Applies gust after aeroelastic convergence
% Michael A. D. Melville, Denver, CO, 80218
if valGUSTTIME > 1
    
    n = 1;
%     valDELTIME = valSDELTIME;        
    [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt,...
        vecLSM, valSPAN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE, matDEF, matTWIST, valNSELE);
    
    matDEF_old = matDEF;
    matTWIST_old = matTWIST;
    
% Runs structure code until static aeroleastic convergence
else
    
    n = 10000;
    
    for tempTIME = 1:n
        
    [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND_STAGGER(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt,...
        vecLSM, valSPAN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE, matDEF, matTWIST, valNSELE, tempTIME);
    
    end
    
    matDEF_old = matDEF;
    matTWIST_old = matTWIST;
%     matDEFGLOB(valTIMESTEP,:) = matDEF(end,:);
%     matTWISTGLOB(valTIMESTEP,:) = matTWIST(end,:);
    
end

[matNPVLST, matNPNEWWAKE, matNEWWAKE, valUINF] = fcnMOVEFLEXWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, vecDVEHVSPN, vecDVELE, matNPVLST, matDEFGLOB,...
    matTWISTGLOB, matSLOPE, valTIMESTEP, vecN, vecM, vecDVEWING, vecDVEPANEL, matSCLST, vecDVEPITCH, matNPDVE, vecSPANDIST, vecCL, valWEIGHT, valAREA, valDENSITY,valUINF);

[ vecDVEHVSPN, vecDVEHVCRD, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEAREA, matDVENORM, matVLST, matDVE, matCENTER, matNEWWAKE ] = fcnVLST2DVEPARAM( matNPDVE, matNPVLST, matNEWWAKE, vecDVETE );

% Determine % relative change between aeroelastic timesteps
tol_def = (100*abs(matDEFGLOB(valTIMESTEP,end)-matDEFGLOB(valTIMESTEP-valSTIFFSTEPS,end))/matDEFGLOB(valTIMESTEP-valSTIFFSTEPS,end));
tol_twist = (100*abs(matTWISTGLOB(valTIMESTEP,end)-matTWISTGLOB(valTIMESTEP-valSTIFFSTEPS,end))/matTWISTGLOB(valTIMESTEP-valSTIFFSTEPS,end));
% 
% % Add in gust velocities to matUINF if convergence tolerance is met
if (tol_def < 25 && tol_twist < 25) || valGUSTTIME > 1
    [matUINF,zvel] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME, zvel, n);
    [matUINF, gust_vel_old] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matCENTER,gust_vel_old);
    valGUSTTIME = valGUSTTIME + 1;
end

end