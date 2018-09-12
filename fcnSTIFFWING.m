function [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, matDEFGLOB, matTWISTGLOB, valUINF, valGUSTTIME, matUINF, flagSTEADY, gust_vel_old, gust_vel_old_vlst, gust_vel_old_ntvlst, gust_vel_old_npvlst] = fcnSTIFFWING(valALPHA, valBETA, valDELTIME, matVLST,...
    matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecN, valTIMESTEP, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF, valGUSTTIME, valGUSTL, valGUSTAMP, flagGUSTMODE, valGUSTSTART, flagSTEADY, matUINF, matUINF_VLST, matUINF_NTVLST, matUINF_NPVLST, gust_vel_old, gust_vel_old_vlst, gust_vel_old_ntvlst, gust_vel_old_npvlst, vecM)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

matCENTER_old = matCENTER;
matNTVLST_old = matNTVLST;
matNPVLST_old = matNPVLST;
matVLST_old = matVLST;

[matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, valUINF] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecCL, valWEIGHT, valAREA, valDENSITY, valTIMESTEP, valUINF, matUINF, matUINF_VLST, matUINF_NTVLST, matUINF_NPVLST, vecN, vecM);

n = 1;

[matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME, n);
[matUINF_vlst] = fcnFLEXUINF(matVLST_old, matVLST, valDELTIME, n);
[matUINF_ntvlst] = fcnFLEXUINF(matNTVLST_old, matNTVLST, valDELTIME, n);
[matUINF_npvlst] = fcnFLEXUINF(matNPVLST_old, matNPVLST, valDELTIME, n);

if valGUSTTIME > 1 || valTIMESTEP == valGUSTSTART
    
    [matUINF, gust_vel_old] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matCENTER,gust_vel_old);
    [matUINF_vlst, gust_vel_old_vlst] = fcnGUSTWING(matUINF_vlst,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matVLST,gust_vel_old_vlst);
    [matUINF_ntvlst, gust_vel_old_ntvlst] = fcnGUSTWING(matUINF_ntvlst,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matNTVLST,gust_vel_old_ntvlst);
    [matUINF_npvlst, gust_vel_old_npvlst] = fcnGUSTWING(matUINF_npvlst,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matNPVLST,gust_vel_old_npvlst);
    valGUSTTIME = valGUSTTIME + 1;
    
end

matDEFGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+1);

matTWISTGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+1);

end