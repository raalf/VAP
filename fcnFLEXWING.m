function [valDELTIME, matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST,...
    vecSPANDIST, matSC, vecMAC, vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF_old, matTWIST_old, matSLOPE,...
    matNPVLST, matNPNEWWAKE, matNEWWAKE, valUINF, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, matVLST, matDVE, matCENTER, matUINF, valGUSTTIME, flagSTEADY] = fcnFLEXWING(vecDVEHVSPN,...
    vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF, vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL,...
    vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecLIFTDIST, vecMOMDIST, valSPAN, valTIMESTEP, matDEFGLOB, matTWISTGLOB,...
    matSLOPE, valALPHA, valBETA, matVLST, matCENTER, matDVE, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF,...
    flagSTATIC, valSDELTIME, valDELTIME, matDEF, matTWIST, valSTIFFSTEPS, valGUSTTIME, valGUSTAMP, valGUSTL, flagGUSTMODE, flagSTEADY)

% valDELTIME = valSDELTIME;

[matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST, matSC, vecMAC] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
    vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW);

matCENTER_old = matCENTER;

% Conditional statement on whether to use constant lift and
% moment (i.e. for static aeroelasticity) or dynamic loads
% that change each timestep (dynamic aeroelasticity)
if valGUSTTIME > 1
    

    valDELTIME = valSDELTIME;
    flagSTEADY = 2;

%     for tempTIME = 1:ceil(valDELTIME/valSDELTIME)
%         
%         [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND_STAGGER(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt,...
%             vecLSM, vecN, valSPAN, vecDVEHVSPN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE, valDELTIME, tempTIME, matDEF, matTWIST);
%         
%     end
    [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt,...
        vecLSM, vecN, valSPAN, vecDVEHVSPN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE);
    
    matDEF_old = matDEF;
    matTWIST_old = matTWIST;
else
    
    for tempTIME = 1:20000
        
        [vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE] = fcnWINGTWISTBEND_STAGGER(vecLIFTDIST, vecMOMDIST, matEIx, vecLM, vecJT, matGJt,...
            vecLSM, vecN, valSPAN, vecDVEHVSPN, valTIMESTEP, matDEFGLOB, matTWISTGLOB, vecSPANDIST, valSDELTIME, matSLOPE, valDELTIME, tempTIME, matDEF, matTWIST);
    end
    matDEF_old = matDEF;
    matTWIST_old = matTWIST;
    matDEFGLOB(valTIMESTEP,:) = matDEF(end,:);
    matTWISTGLOB(valTIMESTEP,:) = matTWIST(end,:);
end

[matNPVLST, matNPNEWWAKE, matNEWWAKE, valUINF] = fcnMOVEFLEXWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, vecDVEHVSPN, vecDVELE, matNPVLST, matDEFGLOB,...
    matTWISTGLOB, matSLOPE, valTIMESTEP, vecN, vecM, vecDVEWING, vecDVEPANEL, matSCLST, vecDVEPITCH, matNPDVE, vecSPANDIST, vecCL, valWEIGHT, valAREA, valDENSITY,valUINF);

[ vecDVEHVSPN, vecDVEHVCRD, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEAREA, matDVENORM, matVLST, matDVE, matCENTER, matNEWWAKE ] = fcnVLST2DVEPARAM( matNPDVE, matNPVLST, matNEWWAKE, vecDVETE );

[matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME, valTIMESTEP);

% Add in gust velocities to matUINF
tol = (100*abs(matDEFGLOB(valTIMESTEP,end-2)-matDEFGLOB(valTIMESTEP-valSTIFFSTEPS,end-2))/matDEFGLOB(valTIMESTEP-valSTIFFSTEPS,end-2));
if tol < 0.05 || valGUSTTIME > 1
    [matUINF] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF);
    valGUSTTIME = valGUSTTIME + 1;
end

end