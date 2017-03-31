function [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, matUINF, matDEFGLOB, matTWISTGLOB, valUINF] = fcnSTIFFWING_STATIC(valALPHA, valBETA,...
    valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecN, valTIMESTEP, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF, matNPDVE,...
    matDEFGLOB, matTWISTGLOB)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

matCENTER_old = matCENTER;

[matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, valUINF] = fcnMOVEWING_STATIC(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST, matNPVLST, vecCL, valWEIGHT, valAREA, valDENSITY, valTIMESTEP, valUINF, matNPDVE);

[matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME, valTIMESTEP);

if valTIMESTEP == 1
    matDEFGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+1);
    matTWISTGLOB(valTIMESTEP,:) = zeros(1,sum(vecN,1)+1);
else
    matDEFGLOB(valTIMESTEP,:) = matDEFGLOB(valTIMESTEP - 1,:);
    matTWISTGLOB(valTIMESTEP,:) = matTWISTGLOB(valTIMESTEP - 1,:);
end

end