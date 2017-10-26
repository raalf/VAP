function [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, matDEFGLOB, matTWISTGLOB, valUINF, valGUSTTIME, matUINF, flagSTEADY, gust_vel_old] = fcnSTIFFWING(valALPHA, valBETA, valDELTIME, matVLST,...
    matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecN, valTIMESTEP, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF, valGUSTTIME, valGUSTL, valGUSTAMP, flagGUSTMODE, valGUSTSTART, flagSTEADY, matUINF, gust_vel_old,test)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

matCENTER_old = matCENTER;

[matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, valUINF] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecCL, valWEIGHT, valAREA, valDENSITY, valTIMESTEP, valUINF, matUINF);

% [matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME);

if valGUSTTIME > 1 || valTIMESTEP == valGUSTSTART
    
    flagSTEADY = 2;
    
    [matUINF, gust_vel_old] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matCENTER,gust_vel_old);
    valGUSTTIME = valGUSTTIME + 1;
    
end

matDEFGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+1);

matTWISTGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+1);

end