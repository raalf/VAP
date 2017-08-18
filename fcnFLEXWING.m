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

[matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST, matSC, vecMAC] = fcnIMPSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
    vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW, matCENTER, valSPAN, vecDVEHVCRD);

matCENTER_old = matCENTER;

% Applies gust after aeroelastic convergence
% Michael A. D. Melville, Denver, CO, 80218
% if valGUSTTIME > 1
    

    valDELTIME = valSDELTIME;
    flagSTEADY = 2;
    valSTIMESTEP = 1;
    valDELTIME = 0.05;

%     for tempTIME = 1:ceil(valDELTIME/valSDELTIME)
%         
%         [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND_STAGGER(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt,...
%             vecLSM, vecN, valSPAN, vecDVEHVSPN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE, valDELTIME, tempTIME, matDEF, matTWIST);
%         
%     end
    [matDEF, matTWIST] = fcnWINGTWISTBEND(valDENSITY,valDELTIME,valSPAN,valAREA,valSTIMESTEP,vecDVEHVSPN,vecDVEHVCRD,...
        vecLEDVES,vecLSAC,vecJT,vecLSM,vecLAMBDA,vecLIFTDIST,vecMOMDIST,valUINF,matEIx,matGJt,matDEF,matTWIST,vecLM,matCENTER);
    
    for valSTIMESTEP = 3:2000
        [matDEF, matTWIST] = fcnWINGTWISTBEND(valDENSITY,valDELTIME,valSPAN,valAREA,valSTIMESTEP,vecDVEHVSPN,vecDVEHVCRD,...
        vecLEDVES,vecLSAC,vecJT,vecLSM,vecLAMBDA,vecLIFTDIST,vecMOMDIST,valUINF,matEIx,matGJt,matDEF,matTWIST,vecLM);
    end
    
    % Use parabolic interpolation to determine deflection and twist at DVE
    % edges rather than control points. This is passed into fcnMOVEFLEXWING
    [matDEFEDGE(valTIMESTEP,:)] = fcnPARABINTERP(matDEF(end,:),matCENTER(vecLEDVES,2),vecSPANDIST(2:end-1));
    [matTWISTEDGE(valTIMESTEP,:)] = fcnPARABINTERP(matDEF(end,:),matCENTER(vecLEDVES,2),vecSPANDIST(2:end-1));
    
    matDEF_old = matDEF;
    matTWIST_old = matTWIST;
    
% Runs structure code until static aeroleastic convergence
% else
    
%     for tempTIME = 1:20000
%         
%         [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND_STAGGER(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt,...
%             vecLSM, vecN, valSPAN, vecDVEHVSPN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE, valDELTIME, tempTIME, matDEF, matTWIST);
%     end
%     matDEF_old = matDEF;
%     matTWIST_old = matTWIST;
%     matDEFGLOB(valTIMESTEP,:) = matDEF(end,:);
%     matTWISTGLOB(valTIMESTEP,:) = matTWIST(end,:);
% end

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