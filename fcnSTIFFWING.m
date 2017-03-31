function [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, matDEFGLOB, matTWISTGLOB, valUINF] = fcnSTIFFWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecN, valTIMESTEP, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

% matCENTER_old = matCENTER;

[matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, valUINF] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecCL, valWEIGHT, valAREA, valDENSITY, valTIMESTEP, valUINF);

% [matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME);

matDEFGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+5);

matTWISTGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+5);

end