function [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, matUINF, matDEFGLOB, matTWISTGLOB] = fcnSTIFFWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecN, valTIMESTEP)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

matCENTER_old = matCENTER;

[matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST);

[matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME);

matDEFGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+5);

matTWISTGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+5);

end