function [valCL,valCLF, valCLI, valCDI, valE, vecDVENFREE, vecDVENIND, vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecLIFTDIST, vecMOMDIST, vecCLDIST, valVINF] = fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, matUINF, vecDVELESWP, vecDVEMCSWP, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE,...
    valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecSYM, vecDVETESWP, valAREA, valSPAN, valBETA, vecDVEWING, vecWDVEWING, vecN, vecM, vecDVEPANEL, vecDVEAREA, vecSPNWSECRD, vecSPNWSEAREA, matQTRCRD, valDENSITY, valWEIGHT, vecLEDVES,...
    vecUINF, matSCLST, vecSPANDIST, matNPVLST, matNPDVE, matSC, vecMAC, valCM, valUINF, matAEROCNTR, flagSTIFFWING,temp)
%% Forces package
%place any force functions in here and add a description.

% INPUT:

% OUTPUT:
% valCL - Lift force coefficient
% valCDI - Ind. drag coefficient
% valE - Span Efficiency

%% Element normal forces, lift forces and side forces (freestream and induced)

[vecDVENFREE, vecDVENIND, vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, matLIFTDIR] = fcnDVENFORCE(matCOEFF...
    ,vecK,matDVE,valNELE,matCENTER,matVLST,matUINF,vecDVELESWP,vecDVEMCSWP,vecDVEHVSPN,vecDVEHVCRD,vecDVEROLL,...
    vecDVEPITCH,vecDVEYAW,vecDVELE,matADJE,valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN,vecWDVEHVCRD,...
    vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecSYM, vecDVETESWP);
%% TE Element induced drag forces

%% Induced Drag force

matUINFTE = matUINF(find(vecDVETE>0),:); % Get trailing edge velocities for de-sweeping in induced drag
[inddrag] =fcnDVEINDDRAG(matCOEFF, matDVE, matVLST, matUINFTE, vecDVEHVSPN, vecDVEHVCRD, vecDVETE, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, ...
    valWSIZE, valTIMESTEP, vecSYM, vecDVEWING, vecWDVEWING, temp);

%% Sum up element forces to generate total wing forces

[valCL, valCLF, valCLI, valCY, valCYF, valCYI, valCDI, valE, vecLIFTDIST, vecMOMDIST, vecCLDIST]= fcnWINGNFORCE(vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, inddrag, vecUINF, valAREA, valSPAN, vecSYM, valBETA, vecDVEAREA,...
    vecDVEWING, vecM, vecN, vecDVELE, vecDVEPANEL, vecSPNWSECRD, vecSPNWSEAREA, matQTRCRD, valDENSITY, valWEIGHT, vecDVEHVSPN, vecLEDVES, matUINF, matSCLST, matCENTER, vecSPANDIST, matNPVLST, matNPDVE, matSC, matLIFTDIR, vecMAC, valCM, valUINF, matAEROCNTR, flagSTIFFWING);
